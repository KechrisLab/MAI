\name{MAI}
\alias{MAI}
\title{Mechanism-Aware Imputation}
\usage{
MAI(data_miss,
    MCAR_algorithm = "BPCA",
    MNAR_algorithm = "Single",
    n_cores = 1,
    assay_ix = 1,
    verbose = TRUE
    )
}
\description{
A two-step approach to imputing missing data in metabolomics.
    Where step 1 predicts missing value mechanisms, and step 2 imputes the
    missing values based on the predicted missing mechanisms.
}
\arguments{
    \item{data_miss}{A matrix or dataframe, or a SummarizedExperiment containing
    missing values designated by "NA" to impute}
    \item{MCAR_algorithm}{The imputation algorithm you wish to use to impute
          MCAR predicted missing values. possible algorithms
          c("BPCA", "Multi_nsKNN", "random_forest")}
    \item{MNAR_algorithm}{The imputation algorithm you wish to use to impute
          MNAR predicted missing values.
          possible algorithms c("Single", "nsKNN")}
    \item{n_cores}{The number of cores you want to utilize. Default is 1 core.
          To use all cores specify n_cores = -1.}
    \item{assay_ix}{If data is a Summarized Experiment then this argument
          defines the index of the assay to impute. Default is set to the first
          assay.}
    \item{verbose}{A toggle to suppress console output. Default is TRUE}
}
\examples{
data(untargeted_LCMS_data)
MAI(data_miss=untargeted_LCMS_data,
    MCAR_algorithm = "BPCA",
    MNAR_algorithm="Single",
    n_cores = 1,
    assay_ix = 1,
    verbose = TRUE)
}
\value{
    When matrix or dataframe returns a list containing the following:
    \item{Imputed Data}{Returns dataframes of MAI imputation as well as
        imputations using the single defined algorithms for MCAR imputation
        and MNAR imputation}
    \item{Estimated Parameters}{Returns the estimated alpha, beta, and gamma
        parameters that define the missingness pattern in the data set}

    When a Summarized Experiment returns:
    \item{Imputed Assay}{Returns the imputed data in the specified assay based
        on the assay_ix assigned}
    \item{Estimated Parameters}{Returns estimated parameters in the metadata of
        the Summarized Experiment as a list}
    \item{Single algorithm imputations}{Returns single algorithm imputations in
        the metadata of the Summarized Experiment as a list}
}
