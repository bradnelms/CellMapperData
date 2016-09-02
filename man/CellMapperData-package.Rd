\name{CellMapperData-package}
\alias{CellMapperData-package}
\alias{CellMapperData}
\docType{package}
\title{
Pre-processed data for use with the CellMapper package
}
\description{
Contains microarray data from several large expression compendia that have been pre-processed for use with the \code{CellMapper} package. All datasets are \code{CellMapperList} objects that were pre-processed using the \code{CMprep} function; these can be provided directly to the \code{CMsearch} function. These pre-processed datasets are recommended for routine searches using the \code{CellMapper} package.
}
\details{
Contains the following six pre-processed datasets:

EH170: pre-processed microarray samples from microdissected human brain regions. Normalized microarray data were downloaded from the Allen Brain Atlas in February 2014 at the link: http://human.brain-map.org/static/download (the 6 files listed under the heading "Complete normalized microarray datasets"). The six normalized datasets were then loaded into R and cancatenated to make a unified expression matrix of 3702 arrays, and pre-processed using the \code{CMprep} function.

EH171: pre-processed microarray samples from the Affymetrix HG_U133PlusV2 platform. The data were obtained from the \code{GSE64985} Bioconductor package and pre-processed using the \code{CMprep} function. The original dataset contained 9395 arrays from diverse biological samples.

EH172: pre-processed microarray samples from the Affymetrix HG_U133A platform. The data were obtained from the \code{E.MTAB.62} Bioconductor package and pre-processed using the \code{CMprep} function. The original dataset contained 5372 arrays from diverse biological samples.

EH173: pre-processed microarray samples from the Affymetrix MG_U74Av2 platform. Normalized microarray data were downloaded from Array Express accession E-MTAB-27 (https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-27/) and processed with the R package \code{bias.0.0.3} to reduce the influence of technical bias (Eklund, et al. 2008). Then Mouse Entrez IDs were then mapped to their corresponding human orthologs as described in Nelms, et al. 2016, and the expression matrix was pre-processed with the \code{CMprep} function. The original dataset contained 1332 arrays from diverse biological samples.

EH174: pre-processed microarray samples from the human intestine. An intestine-specific subset of the 'Engreitz' and 'Lukk' datsets was obtained as described in Section 5.2 of the \code{CellMapper} package vignette. The unprocessed data contained 582 microarrays from the 'Engreitz' dataset and 130 microarrays from the 'Lukk' dataset.

EH175: pre-processed microarray samples from the human kidney. Normalized microarray data were downloaded from the Gene Expression Omnibus from the following accessions: GSE32691, GSE35488, GSE37455, GSE37460, and GSE47185. The five normalized datasets were then loaded into R and cancatenated to make a unified expression matrix of 463 arrays, and pre-processed using the \code{CMprep} function.

}
\author{
Brad Nelms

Maintainer: Brad Nelms <bnelms.research@gmail.com>
}
\references{
B.D. Nelms, L. Waldron, L.A. Barrera, A.W. Weflen, J.A. Goettel, G. Guo, R.K. Montgomery, M.R. Neutra, D.T. Breault, S.B. Snapper, S.H. Orkin, M.L. Bulyk, C. Huttenhower and W.I. Lencer. CellMapper: rapid and accurate inference of gene expression in difficult-to-isolate cell types. Genome Biol 2016, Accepted in Principle.

Hawrylycz MJ, Lein ES, Guillozet-Bongaarts AL, Shen EH, Ng L, Miller J a, van de
Lagemaat LN, Smith K A, Ebbert A, Riley ZL, Abajian C, Beckmann CF, Bernard A, Bertagnolli
D, Boe AF, Cartagena PM, Chakravarty MM, Chapin M, Chong J, Dalley R A, Daly BD, Dang C,
Datta S, Dee N, Dolbeare T a, Faber V, Feng D, Fowler DR, Goldy J, Gregor BW, et al. An
anatomically comprehensive atlas of the adult human brain transcriptome. Nature 2012,
489:391-9.

Engreitz JM, Daigle BJ, Marshall JJ, Altman RB. Independent component analysis: mining microarray data for fundamental human gene expression modules. J Biomed Inform 2010, 43:932-44.

Lukk M, Kapushesky M, Nikkil\"a J, Parkinson H, Goncalves A, Huber W, Ukkonen E, Brazma
A. A global map of human gene expression. Nat Biotechnol 2010, 28:322-4.

Zheng-Bradley X, Rung J, Parkinson H, Brazma A. Large scale comparison of global
gene expression patterns in human and mouse. Genome Biol 2010, 11:R124.

Eklund AC, Szallasi Z. Correction of technical bias in clinical microarray data improves concordance with known biological information. Genome Biol 2008, 9:R26.
}
\seealso{
\code{\link[CellMapper]{CMsearch}}, \code{\link[CellMapper]{CMprep}}
}
\examples{
library(ExperimentHub)
hub <- ExperimentHub()
x <- query(hub, "CellMapperData")
x
\dontrun{
## download Brain Atlas resource (EH170)
BrainAtlas <- x[[1]]
}}
\keyword{datasets, CellMapper}
