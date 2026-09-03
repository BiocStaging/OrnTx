#' Check which RNA-seq tools are available
#' @description Checks fastp, STAR, fasterq-dump availability.
#' @return Named logical vector of tool availability.
#' @export
#' @name checkTools
#' @examples
#' checkTools()
checkTools <- function() {
  tools <- c("fastp", "STAR", "fasterq-dump", "featureCounts")
  available <- vapply(tools, function(t) {
    ret <- tryCatch(
      system2(t, "--version", stdout=TRUE, stderr=TRUE),
      error=function(e) NULL, warning=function(w) NULL)
    !is.null(ret)
  }, logical(1))
  cat("\n=== OrnAtlas Tool Check ===\n")
  for (i in seq_along(available)) {
    status <- if (available[i]) "FOUND" else "NOT FOUND"
    cat(sprintf("  %-20s %s\n", names(available)[i], status))
  }
  if (!any(available)) {
    message("No tools found. Pipeline functions require Linux/macOS/HPC.")
    message("importCustomAtlas() and all analysis functions work on all platforms.")
  }
  invisible(available)
}


#' Download RNA-seq data from NCBI SRA
#' @description Downloads SRA files and converts to FASTQ using
#'   fasterq-dump or direct AWS S3 download.
#' @param accessions Character vector of SRA run accessions.
#' @param output_dir Directory to save FASTQ files.
#' @param threads Threads for fasterq-dump (default 4).
#' @param compress Compress output with gzip (default TRUE).
#' @param method "fasterq-dump" or "aws" (direct S3, faster).
#' @return Invisibly returns paths to downloaded FASTQ files.
#' @export
#' @name downloadSRA
#' @examples
#' \dontrun{
#' downloadSRA(
#'   accessions = c("SRR10160942", "SRR10160944"),
#'   output_dir = "raw_reads/",
#'   threads    = 4
#' )
#' }
downloadSRA <- function(accessions, output_dir="raw_reads",
                         threads=4L, compress=TRUE,
                         method=c("aws","fasterq-dump")) {
  method <- match.arg(method)
  dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
  message("[OrnAtlas] Downloading ", length(accessions), " accessions...")
  downloaded <- character(0)
  for (srr in accessions) {
    message("[OrnAtlas] Downloading: ", srr)
    if (method == "aws") {
      url <- paste0(
        "https://sra-pub-run-odp.s3.amazonaws.com/sra/", srr, "/", srr)
      sra_file <- file.path(output_dir, paste0(srr, ".sra"))
      tryCatch({
        utils::download.file(url, destfile=sra_file, mode="wb")
        args <- c("--split-files", "--outdir", output_dir,
                  "--threads", as.character(threads), sra_file)
        ret <- system2("fasterq-dump", args)
        if (ret == 0 && compress) {
          fq <- list.files(output_dir,
                           pattern=paste0(srr, ".*[.]fastq$"),
                           full.names=TRUE)
          for (f in fq) system2("gzip", f)
          downloaded <- c(downloaded, paste0(fq, ".gz"))
        }
        message("[OrnAtlas] ", srr, " done!")
      }, error=function(e) warning("Failed: ", srr, " - ", e$message))
    } else {
      args <- c("--split-files", "--outdir", output_dir,
                "--threads", as.character(threads), srr)
      ret <- system2("fasterq-dump", args)
      if (ret == 0) {
        if (compress) {
          fq <- list.files(output_dir,
                           pattern=paste0(srr, ".*[.]fastq$"),
                           full.names=TRUE)
          for (f in fq) system2("gzip", f)
        }
        downloaded <- c(downloaded,
                         list.files(output_dir,
                                    pattern=paste0(srr, ".*[.]fastq"),
                                    full.names=TRUE))
        message("[OrnAtlas] ", srr, " done!")
      }
    }
  }
  message("[OrnAtlas] Downloaded ", length(downloaded), " files")
  invisible(downloaded)
}


#' Quality control of FASTQ files using fastp
#' @description Runs fastp QC on FASTQ files. Supports paired and
#'   single-end data. Requires fastp to be installed.
#' @param input_dir Directory containing raw FASTQ files.
#' @param output_dir Directory for clean FASTQ output.
#' @param paired Logical. TRUE for paired-end (default).
#' @param threads Threads per sample (default 4).
#' @param min_length Minimum read length after trimming (default 36).
#' @param quality Minimum base quality score (default 20).
#' @return Invisibly returns paths to cleaned FASTQ files.
#' @export
#' @name runQC
#' @examples
#' \dontrun{
#' runQC(input_dir="raw_reads/", output_dir="clean_reads/")
#' }
runQC <- function(input_dir, output_dir="clean_reads",
                   paired=TRUE, threads=4L,
                   min_length=36L, quality=20L) {
  dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
  dir.create(file.path(output_dir, "qc_reports"),
             showWarnings=FALSE, recursive=TRUE)
  if (paired) {
    r1_files <- list.files(input_dir,
      pattern="_1[.]fastq[.]gz$|_R1[.]fastq[.]gz$", full.names=TRUE)
    samples <- gsub("_1[.]fastq[.]gz$|_R1[.]fastq[.]gz$", "",
                     basename(r1_files))
  } else {
    r1_files <- list.files(input_dir,
      pattern="[.]fastq[.]gz$", full.names=TRUE)
    samples <- gsub("[.]fastq[.]gz$", "", basename(r1_files))
  }
  message("[OrnAtlas] Running QC on ", length(samples), " samples...")
  cleaned <- character(0)
  for (i in seq_along(samples)) {
    s    <- samples[i]
    r1   <- r1_files[i]
    out1 <- file.path(output_dir, paste0(s, "_1_clean.fastq.gz"))
    message("[OrnAtlas] QC: ", s)
    if (paired) {
      r2   <- gsub("_1[.]fastq[.]gz$|_R1[.]fastq[.]gz$",
                   "_2.fastq.gz", r1)
      out2 <- file.path(output_dir, paste0(s, "_2_clean.fastq.gz"))
      args <- c("--in1", r1, "--in2", r2,
                "--out1", out1, "--out2", out2,
                "--json", file.path(output_dir, "qc_reports",
                                    paste0(s, ".json")),
                "--html", file.path(output_dir, "qc_reports",
                                    paste0(s, ".html")),
                "--thread", threads,
                "--qualified_quality_phred", quality,
                "--length_required", min_length, "-y")
    } else {
      args <- c("--in1", r1, "--out1", out1,
                "--json", file.path(output_dir, "qc_reports",
                                    paste0(s, ".json")),
                "--thread", threads,
                "--qualified_quality_phred", quality,
                "--length_required", min_length, "-y")
    }
    ret <- system2("fastp", args)
    if (ret == 0) {
      message("[OrnAtlas] ", s, " QC complete!")
      cleaned <- c(cleaned, out1)
    } else warning("[OrnAtlas] fastp failed for: ", s)
  }
  message("[OrnAtlas] QC complete: ", length(cleaned), " samples")
  invisible(cleaned)
}


#' Build STAR genome index
#' @description Builds a STAR genome index from FASTA and GTF.
#'   Required once per genome. Requires STAR to be installed.
#' @param genome_fasta Path to reference genome FASTA file.
#' @param gtf_file Path to genome annotation GTF file.
#' @param index_dir Output directory for genome index.
#' @param threads Threads for index building (default 16).
#' @param ram_gb Available RAM in GB (default 60).
#' @param sa_index_nbases SA index parameter (default 13).
#' @return Invisibly returns the index directory path.
#' @export
#' @name buildGenomeIndex
#' @examples
#' \dontrun{
#' buildGenomeIndex(
#'   genome_fasta = "rosa_genome.fna",
#'   gtf_file     = "rosa_annotation.gtf",
#'   index_dir    = "star_index/"
#' )
#' }
buildGenomeIndex <- function(genome_fasta, gtf_file,
                              index_dir="star_index",
                              threads=16L, ram_gb=60,
                              sa_index_nbases=13L) {
  if (!file.exists(genome_fasta))
    stop("Genome FASTA not found: ", genome_fasta)
  if (!file.exists(gtf_file))
    stop("GTF not found: ", gtf_file)
  dir.create(index_dir, showWarnings=FALSE, recursive=TRUE)
  message("[OrnAtlas] Building STAR genome index...")
  message("[OrnAtlas] Genome: ", basename(genome_fasta))
  message("[OrnAtlas] This may take 20-40 minutes...")
  args <- c("--runMode", "genomeGenerate",
            "--genomeDir", index_dir,
            "--genomeFastaFiles", genome_fasta,
            "--sjdbGTFfile", gtf_file,
            "--runThreadN", as.character(threads),
            "--genomeSAindexNbases", as.character(sa_index_nbases),
            "--limitGenomeGenerateRAM",
            as.character(as.integer(ram_gb * 1e9)))
  ret <- system2("STAR", args)
  if (ret == 0) {
    message("[OrnAtlas] Genome index built successfully!")
  } else {
    stop("[OrnAtlas] STAR indexing failed. Check installation and RAM.")
  }
  invisible(index_dir)
}


#' Align RNA-seq reads using STAR
#' @description Aligns FASTQ files to a reference genome using STAR
#'   and generates ReadsPerGene.out.tab count files compatible with
#'   importCustomAtlas(). Requires STAR to be installed.
#' @param fastq_dir Directory containing cleaned FASTQ files.
#' @param index_dir Path to STAR genome index directory.
#' @param output_dir Directory for alignment output.
#' @param paired Logical. TRUE for paired-end (default).
#' @param threads Threads per sample (default 16).
#' @param ram_gb RAM for BAM sorting in GB (default 60).
#' @param strandedness Library strandedness (default "reverse").
#' @return Invisibly returns paths to count files.
#' @export
#' @name alignReads
#' @examples
#' \dontrun{
#' alignReads(
#'   fastq_dir  = "clean_reads/",
#'   index_dir  = "star_index/",
#'   output_dir = "aligned/"
#' )
#' }
alignReads <- function(fastq_dir, index_dir,
                        output_dir="aligned",
                        paired=TRUE, threads=16L,
                        ram_gb=60, strandedness="reverse") {
  dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
  if (paired) {
    r1_files <- list.files(fastq_dir,
      pattern="_1_clean[.]fastq[.]gz$|_R1_clean[.]fastq[.]gz$",
      full.names=TRUE)
    samples <- gsub("_1_clean[.]fastq[.]gz$|_R1_clean[.]fastq[.]gz$",
                     "", basename(r1_files))
  } else {
    r1_files <- list.files(fastq_dir,
      pattern="_clean[.]fastq[.]gz$", full.names=TRUE)
    samples <- gsub("_clean[.]fastq[.]gz$", "", basename(r1_files))
  }
  message("[OrnAtlas] Aligning ", length(samples), " samples...")
  count_files <- character(0)
  for (i in seq_along(samples)) {
    s    <- samples[i]
    r1   <- r1_files[i]
    sdir <- file.path(output_dir, s)
    dir.create(sdir, showWarnings=FALSE)
    count_file <- file.path(sdir, "ReadsPerGene.out.tab")
    if (file.exists(count_file)) {
      message("[OrnAtlas] ", s, " already aligned, skipping.")
      count_files <- c(count_files, count_file)
      next
    }
    message("[OrnAtlas] Aligning: ", s)
    reads_in <- if (paired) {
      r2 <- gsub("_1_clean[.]fastq[.]gz$|_R1_clean[.]fastq[.]gz$",
                 "_2_clean.fastq.gz", r1)
      c(r1, r2)
    } else r1
    args <- c("--runThreadN", as.character(threads),
              "--genomeDir", index_dir,
              "--readFilesIn", reads_in,
              "--readFilesCommand", "zcat",
              "--outSAMtype", "BAM", "SortedByCoordinate",
              "--quantMode", "GeneCounts",
              "--outFileNamePrefix", paste0(sdir, "/"),
              "--limitBAMsortRAM",
              as.character(as.integer(ram_gb * 1e9)))
    ret <- system2("STAR", args)
    if (ret == 0) {
      message("[OrnAtlas] ", s, " aligned!")
      count_files <- c(count_files, count_file)
    } else warning("[OrnAtlas] STAR failed for: ", s)
  }
  message("[OrnAtlas] Alignment complete: ",
          length(count_files), " count files")
  invisible(count_files)
}


#' Run complete RNA-seq pipeline from FASTQ to SummarizedExperiment
#' @description Convenience wrapper running the complete OrnAtlas
#'   pipeline: QC (fastp) -> Alignment (STAR) -> Import -> Normalize.
#'   For SRA download first, use downloadSRA() separately.
#'   Requires fastp and STAR installed (Linux/macOS/HPC only).
#' @param fastq_dir Directory containing raw FASTQ files.
#' @param genome_fasta Path to reference genome FASTA.
#' @param gtf_file Path to genome annotation GTF.
#' @param metadata_file Path to sample metadata CSV/TSV.
#' @param species Scientific name of the plant species.
#' @param work_dir Working directory for intermediate files.
#' @param threads Threads for alignment (default 16).
#' @param strandedness Library strandedness (default "reverse").
#' @param paired Logical. Paired-end data (default TRUE).
#' @param skip_qc Logical. Skip fastp QC step (default FALSE).
#' @param skip_index Logical. Skip index building if built.
#' @return A CPM-normalized SummarizedExperiment.
#' @export
#' @name runFullPipeline
#' @examples
#' \dontrun{
#' se <- runFullPipeline(
#'   fastq_dir     = "raw_reads/",
#'   genome_fasta  = "rosa_genome.fna",
#'   gtf_file      = "rosa_annotation.gtf",
#'   metadata_file = "metadata.csv",
#'   species       = "Rosa chinensis",
#'   threads       = 16
#' )
#' plotPCA(se, color_by = "tissue")
#' }
runFullPipeline <- function(fastq_dir,
                             genome_fasta,
                             gtf_file,
                             metadata_file,
                             species       = "Unknown species",
                             work_dir      = "ornAtlas_pipeline",
                             threads       = 16L,
                             strandedness  = "reverse",
                             paired        = TRUE,
                             skip_qc       = FALSE,
                             skip_index    = FALSE) {
  message("\n========================================")
  message("  OrnAtlas Full Pipeline")
  message("  Species: ", species)
  message("========================================\n")
  tools <- checkTools()
  if (!tools["STAR"])
    stop("STAR not found. Install STAR and add to PATH.")
  if (!skip_qc && !tools["fastp"])
    stop("fastp not found. Install fastp or set skip_qc=TRUE.")
  clean_dir  <- file.path(work_dir, "clean_reads")
  index_dir  <- file.path(work_dir, "star_index")
  align_dir  <- file.path(work_dir, "aligned")
  dir.create(work_dir, showWarnings=FALSE, recursive=TRUE)
  if (!skip_qc) {
    message("\n--- Step 1: Quality Control ---")
    runQC(input_dir=fastq_dir, output_dir=clean_dir,
          paired=paired, threads=min(threads, 4L))
    fastq_input <- clean_dir
  } else {
    message("\n--- Step 1: Skipping QC ---")
    fastq_input <- fastq_dir
  }
  if (!skip_index) {
    message("\n--- Step 2: Building STAR Index ---")
    buildGenomeIndex(genome_fasta=genome_fasta,
                     gtf_file=gtf_file,
                     index_dir=index_dir,
                     threads=threads)
  } else {
    message("\n--- Step 2: Using existing index ---")
  }
  message("\n--- Step 3: STAR Alignment ---")
  alignReads(fastq_dir=fastq_input, index_dir=index_dir,
             output_dir=align_dir, paired=paired,
             threads=threads, strandedness=strandedness)
  message("\n--- Step 4: Importing counts ---")
  se <- importCustomAtlas(
    counts_dir    = align_dir,
    metadata_file = metadata_file,
    species       = species,
    strandedness  = strandedness
  )
  message("\n--- Step 5: Normalizing ---")
  se <- normalizeCounts(se, method="CPM")
  message("\n========================================")
  message("  Pipeline Complete!")
  message("  Genes:   ", nrow(se))
  message("  Samples: ", ncol(se))
  message("  Ready for plotPCA(), runDE(), getFlowerGenes()")
  message("========================================\n")
  se
}
