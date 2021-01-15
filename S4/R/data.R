#' Microarray data of leukemia after preprocessing(GSE6891)
#'
#' Expression matrix with 89 samples(Row) and 2000 genes(columns). Originally there are 54,676 probesets in each dataset
#' (GSE6891, GSE17855 and GSE13159) and we remove the probesets with missing values and select the probesets with the largest interquartile
#' range to represent the gene if multiple probesets are mapped to the same gene. 20,192 unique genes
#' remained for every study after this preprocessing. Furthermore, for each study, we transform data
#' to log scale and only keep the top 10,000 genes with the largest mean expression level (i.e. filter out low-expressed genes).
#' We further lter out 8,000 genes with smaller variance (i.e. genes with little predictive information).
#' Finally, the remaining 2,000 genes are used in the analysis.
#' @format a list of two components. Expression matrix with 89 samples(Row) and 2000 genes(columns); A vector of class label for samples.
#' @references Verhaak, R. G., Wouters, B. J., Erpelinck, C. A., Abbas, S., Beverloo, H. B., Lugthart, S., ... and Valk, P. J. (2009).
#' Prediction of molecular subtypes in acute myeloid leukemia based on gene expression profiling. haematologica, 94(1), 131-134.
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6891}
"ds.GSE6891"



#' Microarray data of leukemia after preprocessing(GSE17855)
#'
#' Expression matrix with 74 samples(Row) and 2000 genes(columns). Originally there are 54,676 probesets in each dataset
#' (GSE6891, GSE17855 and GSE13159) and we remove the probesets with missing values and select the probesets with the largest interquartile
#' range to represent the gene if multiple probesets are mapped to the same gene. 20,192 unique genes
#' remained for every study after this preprocessing. Furthermore, for each study, we transform data
#' to log scale and only keep the top 10,000 genes with the largest mean expression level (i.e. filter out low-expressed genes).
#' We further lter out 8,000 genes with smaller variance (i.e. genes with little predictive information).
#' Finally, the remaining 2,000 genes are used in the analysis.
#' @format a list of two components. Expression matrix with 74 samples(Row) and 2000 genes(columns); A vector of class label for samples.
#' @references Balgobind, B. V., Van den Heuvel-Eibrink, M. M., Menezes, R. X., Reinhardt, D., Hollink, I. H.,
#' Peters, S. T., van Wering, E. R., Kaspers, G. J., Cloos, J., de Bont, E. S. M., et al. (2010).
#' Evaluation of gene expression signatures predictive for cytogenetic and molecular subtypes of
#' pediatric acute myeloid leukemia. haematologica, pages haematol-2010.
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17855}
"ds.GSE17855"


#' Microarray data of Leukemia after preprocessing(GSE13159)
#'
#' Expression matrix with 105 samples(Row) and 2000 genes(columns). Originally there are 54,676 probesets in each dataset
#' (GSE6891, GSE17855 and GSE13159) and we remove the probesets with missing values and select the probesets with the largest interquartile
#' range to represent the gene if multiple probesets are mapped to the same gene. 20,192 unique genes
#' remained for every study after this preprocessing. Furthermore, for each study, we transform data
#' to log scale and only keep the top 10,000 genes with the largest mean expression level (i.e. filter out low-expressed genes).
#' We further lter out 8,000 genes with smaller variance (i.e. genes with little predictive information).
#' Finally, the remaining 2,000 genes are used in the analysis.
#' @format A list of two components: Expression matrix with 105 samples(Row) and 2000 genes(columns); a vector of class label for samples.
#' @references Kohlmann, A., Kipps, T. J., Rassenti, L. Z., Downing, J. R., Shurtleff, S. A., Mills, K. I., Gilkes,
#' A. F., Hofmann, W.-K., Basso, G., DellOrto, M. C., et al. (2008). An international standardization
#' programme towards the application of gene expression profiling in routine leukaemia diagnostics:
#'  the microarray innovations in leukemia study prephase. British journal of haematology,142(5):802-807.
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13159}
"ds.GSE13159"

#' Microarray data after preprocessing(Mammalian tissue types dataset:)
#'
#' Expression matrix with 102 samples(Row) and 2000 genes(columns). Gene expression from human and mouse samples across a diverse
#' array of tissues, organs, and cell lines have been profiled by Su et al. (2002). Here we
#' only consider four tissue types: breast, prostate, lung, and colon.
#' The original data has 102 samples and 5565 probesets (genes). Following
#' similar proprocessing procedure, We filter 3000 genes with highest expression value and then 2000
#' genes were used in the analysis after further filtering low-variance genes.
#' @format A list of two components: Expression matrix with 102 samples(Row) and 2000 genes(columns); a vector of class label for samples.
#' @references Su, A. I., Cooke, M. P., Ching, K. A., Hakak, Y., Walker, J. R., Wiltshire, T., Orth, A. P., Vega,
#' R. G., Sapinoso, L. M., Moqrich, A., et al. (2002). Large-scale analysis of the human and mouse
#' transcriptomes. Proceedings of the National Academy of Sciences, 99(7):4465-4470.
#' @source \url{http://portals.broadinstitute.org/cgi-bin/cancer/datasets.cgi}
"ds.TissueType"

#' RNA Sequencing data of rat brain after preprocessing(GSE47474)
#'
#' Expression matrix with 36 samples(Row) and 2000 genes(columns). RNA samples from three brain regions (hippocampus, striatum and
#' prefrontal cortex) were sequenced for both control strains and HIV strains. Only the 36 control
#' strains (12 in each brain region) are used here to see whether samples from three brain regions can
#' be correctly clustered The original count data is transformed into
#' CPM value and then 2000 genesare kept by filtering low-expressed genes and low-variance genes.
#' @format A list of two components: Expression matrix with 36 samples(Row) and 2000 genes(columns); a vector of class label for samples.
#' @references Li, M. D., Cao, J., Wang, S., Wang, J., Sarkar, S., Vigorito, M., Ma, J. Z., and Chang, S. L. (2013).
#' Transcriptome sequencing of gene expression in the brain of the hiv-1 transgenic rat. PLoS One, 8(3):e59582.
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47474}
"ds.GSE47474"

#' RNA Sequencing data after preprocessing(Pancancer)
#'
#' Expression matrix with 801 samples(Row) and 2000 genes(columns). This collection of data consists of five different
#' types of tumor: 300 breast cancer (BRCA), 146 kidney clear cell carcinoma (KIRC), 78 colon
#' cancer (COAD), 141 lung adenocarcinoma (LUAD) and 136 prostate cancer (PRAD). The data
#' has already been normalized and we use the same filtering process to keep 2000 genes; a vector of class label for samples.
#' @format A list of two components: Expression matrix with 801 samples(Row) and 2000 genes(columns).
#' @source \url{https://archive.ics.uci.edu/ml/datasets/gene+expression+cancer+RNA-Seq}
"ds.Pancancer"

#' SNP data after preprocessing
#'
#' Expression matrix with 293 samples(Row) and 17026 SNP(columns). Following the same preprocessing procedure as Witten
#' and Tibshirani (2010), only phase III SNP data is used and we restrict the analysis to chromosome
#' 22 of three populations: African ancestry in southwest USA (ASW), Utah residents with European
#' ancestry (CEU), and Han Chinese from Beijing (CHB) since these three populations are known 26
#' to be genetically distinct. All the available SNPs on chromesome22 are considered in the data,
#' which gives us 293 samples and 17026 SNP. We then coded AA as 2, Aa as 1 and aa as 0, and use
#' 5-nearest neighbors method(Troyanskaya et al., 2001) to impute the missing data.
#' @format A list of two components: Matrix with 293 samples(Row) and 17026 genes(columns); a vector of class label for samples.
#'
#' @references (1): Troyanskaya, O., Cantor, M., Sherlock, G., Brown, P., Hastie, T., Tibshirani, R., Botstein, D., and
#' Altman, R. B. (2001). Missing value estimation methods for dna microarrays. Bioinformatics, 17(6):520-525. (2):
#' Witten, D. M. and Tibshirani, R. (2010). A framework for feature selection in clustering. Journal
#  of the American Statistical Association, 105(490):713-726.

#' @source \url{ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2008-07_phaseIII/hapmap_format/forward/}
"ds.SNP"


#' Plant species leaves dataset after preprocessing
#'
#' Matrix with 64 samples(Row) and 190 features(columns). Mallah et al. (2013) introduced a dataset consisting of one-hundred
#'species of plants with three types of features for leaves: shape texture and margin. Here we
#'only consider 4 species out of 100, alnus rubra, acer capillipes, cornus macrophylla and quercus
#'chrysolepis. After deleting features with any missing values, we have 64 samples (16 for each species) and 190 features.
#' @format A list of two components: Matrix with 64 samples(Row) and 190 features(columns); a vector of class label for samples.
#' @references Mallah, C., Cope, J., and Orwell, J. (2013). Plant leaf classification using probabilistic integration
#' of shape, texture and margin features. Signal Processing, Pattern Recognition and Applications, 5(1).
#' @source \url{https://archive.ics.uci.edu/ml/datasets/One-hundred+plant+species+leaves+data+set}
"ds.Leaf"

#' letter-name dataset after preprocessing
#'
#' Matrix with 617 samples(Row) and 1200 features(columns). ISOLET dataset was generated by a study where 150 subjects spoke each letter
#' of the alphabet twice and recorded 617 features including spectral coefficients, contour features,
#' sonorant features, pre-sonorant features and post-sonorant features. We only use five vowels and
#' 1200 training subjects ( 240 samples for each of five vowels).
#' @format A list of two components: Matrix with 617 samples(Row) and 1200 features(columns); a vector of class label for samples.
#' @references Dheeru, D. and Karra Taniskidou, E. (2017). UCI machine learning repository.
#' @source \url{https://archive.ics.uci.edu/ml/datasets/isolet}
"ds.ISOLET"
