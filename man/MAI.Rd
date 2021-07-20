\name{MAI}
\alias{MAI}
\title{Mechanism-Aware Imputation}
\usage{
MAI(data_miss,
    MCAR_algorithm = "BPCA",
    MNAR_algorithm = "Single",
    n_cores = 1
    )
}
\description{
A two-step approach to imputing missing data in metabolomics.
    Where step 1 predicts missing value mechanisms, and step 2 imputes the
    missing values based on the predicted missing mechanisms.
}
\arguments{
    \item{data_miss}{A matrix containing missing values to impute}
    \item{MCAR_algorithm}{The imputation algorithm you wish to use to impute
          MCAR predicted missing values.}
    \item{MNAR_algorithm}{The imputation algorithm you wish to use to impute
          MNAR predicted missing values.}
    \item{n_cores}{The number of cores you want to utilize. Default is 1 core.
          To use all cores specify n_cores = -1.}
}
\examples{
data(data_miss)
MAI(data_miss=data_miss,
    MCAR_algorithm = "BPCA",
    MNAR_algorithm="Single",
    n_cores = 1)
}
\value{
    \item{Imputed Data}{Returns dataframes of MAI imputation as well as
        imputations using the single defined algorithms for MCAR imputation
        and MNAR imputation}
    \item{Estimated Parameters}{Returns the estimated alpha, beta, and gamma
        parameters that define the missingness pattern in the data set}
}
