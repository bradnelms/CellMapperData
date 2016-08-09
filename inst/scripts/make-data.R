##########################################
## EH170 Dataset
##########################################

# Normalized microarray data were downloaded from the Allen Brain Atlas in February 2014 at the link: http://human.brain-map.org/static/download (the 6 files listed under the heading "Complete normalized microarray datasets"). The six normalized datasets were then loaded into R and cancatenated to make a unified expression matrix of 3702 arrays, and pre-processed using the \code{CM.prep} function.



##########################################
## EH171 Dataset
##########################################

library(ExperimentHub)
library(CellMapper)
hub <- ExperimentHub()
x <- query(hub, "HumanAffyData")
x <- x[[1]]
EH171 <- CM.prep(x, DataSource = "Engreitz, et al. (2010)", GeneIDType = "Human Entrez IDs")



##########################################
## EH172 Dataset
##########################################

library(ExperimentHub)
library(CellMapper)
hub <- ExperimentHub()
x <- query(hub, "HumanAffyData")
x <- x[[2]]
EH172 <- CM.prep(x, DataSource = "Lukk, et al. (2010)", GeneIDType = "Human Entrez IDs")



##########################################
## EH173 Dataset
##########################################

# Normalized microarray data were downloaded from Array Express accession E-MTAB-27 (https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-27/) and processed with the R package \code{bias.0.0.3} to reduce the influence of technical bias (Eklund, et al. 2008). Then Mouse Entrez IDs were then mapped to their corresponding human orthologs as described in Nelms, et al. 2016, and the expression matrix was pre-processed with the \code{CM.prep} function.



##########################################
## EH174 Dataset
##########################################

library(ExperimentHub)
library(CellMapper)
hub <- ExperimentHub()
x <- query(hub, "HumanAffyData")

EH174 <- list()
E.MTAB.62 <- x[[2]]
pDat <- pData(E.MTAB.62)
samples <- which(pDat$OrganismPart %in% c("colon", "colon mucosa", "Small intestine"))
Lukk_unprocessed.gut <- exprs(E.MTAB.62)[, samples] 
EH174[['Lukk.Gut']] <- CM.prep(Lukk_unprocessed.gut, DataSource = "Gut-specific subset of samples from Lukk, et al. (2010)", GeneIDType = "Human Entrez IDs")

GSE64985 <- x[[1]]
pDat <- pData(GSE64985)
keywords <- c("colon", "intestin")
select <- grepl(paste(keywords, collapse = "|"), pDat$title, ignore.case = TRUE) | grepl(paste(keywords, collapse = "|"), pDat$description, ignore.case = TRUE)
Engreitz_unprocessed.gut <- exprs(GSE64985)[,select]
EH174[['Engreitz.gut']] <- CM.prep(Engreitz_unprocessed.gut, DataSource = "Gut-specific subset of samples from Engreitz, et al. (2010)", GeneIDType = "Human Entrez IDs")



##########################################
## EH175 Dataset
##########################################

# Normalized microarray data were downloaded from the Gene Expression Omnibus from the following accessions: GSE32691, GSE35488, GSE37455, GSE37460, and GSE47185. The five normalized datasets were then loaded into R and cancatenated to make a unified expression matrix of 463 arrays, and pre-processed using the \code{CM.prep} function.

