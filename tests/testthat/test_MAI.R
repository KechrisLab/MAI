library(MAI)
data("untargeted_LCMS_data")

test_that("MAI works properly", {
    result = MAI(untargeted_LCMS_data)
    expect_true(ncol(result[["Imputed_data"]])==ncol(untargeted_LCMS_data))
    expect_true(sum(is.na(result[["Imputed_data"]]))==0)
    expect_true(is.numeric(result[["Imputed_data"]]))
    expect_true(length(result[["Estimated_Params"]])==3)
    expect_true(is.numeric(result[["Estimated_Params"]][["Alpha"]]))
    expect_true(is.numeric(result[["Estimated_Params"]][["Beta"]]))
    expect_true(is.numeric(result[["Estimated_Params"]][["Gamma"]]))
})

