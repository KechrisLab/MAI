library(MAI)
data("data_miss")

test_that("MAI works properly", {
    result = MAI(data_miss)
    expect_true(ncol(result[["Imputed_data"]][["MAI"]])==ncol(data_miss))
    expect_true(sum(is.na(result[["Imputed_data"]][["MAI"]]))==0)
    expect_true(is.numeric(result[["Imputed_data"]][["MAI"]]))
    expect_true(length(result[["Estimated_Params"]])==3)
    expect_true(is.numeric(result[["Estimated_Params"]][["Alpha"]]))
    expect_true(is.numeric(result[["Estimated_Params"]][["Beta"]]))
    expect_true(is.numeric(result[["Estimated_Params"]][["Gamma"]]))
})

