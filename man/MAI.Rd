\name{MAI}
\alias{MAI}
\title{Mechanism-Aware Imputation}
\usage{
MAI(data_miss,
    MCAR_algorithm = c("BPCA", "Multi_nsKNN", "random_forest"),
    MNAR_algorithm = c("nsKNN", "Single"),
    n_cores = 1,
    assay_ix = 1,
    forest_list_args = list(
        ntree = 300,
        proximity = FALSE
    ),
    verbose = TRUE
    )
}
\description{
A two-step approach to imputing missing data in metabolomics.
    Step 1 uses a random forest classifier to classify missing values as
    either Missing Completely at Random/Missing At Random (MCAR/MAR) or Missing
    Not At Random (MNAR). MCAR/MAR are combined because it is often difficult to
    distinguish these two missing types in metabolomics data. Step 2 imputes the
    missing values based on the classified missing mechanisms, using the
    appropriate imputation algorithms. Imputation algorithms tested and
    available for MCAR/MAR include Bayesian Principal Component Analysis (BPCA),
    Multiple Imputation No-Skip K-Nearest Neighbors (Multi_nsKNN), and
    Random Forest. Imputation algorithms tested and available for MNAR include
    nsKNN and a single imputation approach for imputation of metabolites where
    left-censoring is present.
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
    \item{forest_list_args}{Random forest named arguments to pass to the random forest training process. Defualt args are ntree = 300 and proximity = FALSE}
    \item{verbose}{A toggle to suppress console output. Default is TRUE}
}
\examples{
data(untargeted_LCMS_data)
MAI(data_miss=untargeted_LCMS_data,
    MCAR_algorithm = "BPCA",
    MNAR_algorithm="Single",
    n_cores = 1,
    assay_ix = 1,
    forest_list_args = list(
        ntree = 300,
        proximity = FALSE
    ),
    verbose = TRUE)
}
\value{
    When matrix or dataframe returns a list containing the following:
    \item{Imputed Data}{Returns dataframes of MAI imputation}
    \item{Estimated Parameters}{Returns the estimated \eqn{\alpha}, \eqn{\beta},
    and \eqn{\gamma} parameters that define the missingness pattern in the data
    set. These parameters estimate the ratio of MCAR/MAR to MNAR in the data.
    The parameters \eqn{\alpha} and \eqn{\beta} separate high, medium, and low
    average abundance metabolites, while the parameter \eqn{\gamma} is used to
    impose missingness in the medium and low abundance metabolites.  A smaller
    \eqn{\alpha} corresponds to more MCAR/MAR being present, while larger
    \eqn{\beta} and \eqn{\gamma} values imply more MNAR values being present.
    The returned estimated parameters are then used to impose known missingness
    in the complete subset of the input data. Subsequently, a random forest
    classifier is trained to classify the known missingness in the complete
    subset of the input data. Once the classifier is established it is applied
    to the unknown missingness of the full input data to predict the
    missingness. Finally, the missing values are imputed using a specific
    algorithm, chosen by the user, according to the predicted
    missingness mechanism.}

    When a Summarized Experiment returns:
    \item{Imputed Assay}{Returns the imputed data in the specified assay based
        on the assay_ix assigned}
    \item{Estimated Parameters}{Returns estimated parameters in the metadata of
        the Summarized Experiment as a list}
}
