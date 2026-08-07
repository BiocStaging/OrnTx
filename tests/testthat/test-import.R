test_that("importRosaCounts gives helpful error for missing file", {
  expect_error(
    importRosaCounts(
      file_paths = c(Sample1 = "nonexistent.tab"),
      col_data   = data.frame(stage = "bud",
                              row.names = "Sample1")
    ),
    regexp = "cannot open"
  )
})

test_that("buildAtlas requires valid group_var", {
  expect_error(
    buildAtlas(
      SummarizedExperiment::SummarizedExperiment(),
      group_var = "nonexistent"
    )
  )
})
