% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Grid_Enrichment.R
\name{Enrich_Networks}
\alias{Enrich_Networks}
\title{Enrichment analysis for genes in the network modules}
\usage{
Enrich_Networks(gene_name, background, db)
}
\arguments{
\item{gene_name}{Gene symbol for the module, in a form of character vector}

\item{background}{Background gene set user specified}

\item{db}{Database of interest: KEGG or GO}
}
\description{
The function enrichment() is dealing with the multimiR target genes, enrich network here is dealing with
network genes for each network module
}
