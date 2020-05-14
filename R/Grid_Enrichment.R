#' Obtain CV error for all grids
#'
#' @param Dir The working directory from previous running
#' @param grids Number of grids that are using from the previous running
#' @param subSamp The subsampling number used in previous running, must be the same since this is the
#' call that read in the data file
#' @example my_error <- get_error(CVDirs,grids = 20, subSamp = 10)
#' @export

get_error <- function(Dir,grids, subSamp)
{
  error <- matrix(,, ncol = 7)
  sd_train <- rep(NA,grids)
  sd_test <- rep(NA,grids)
  sd_error <- rep(NA,grids)
  for (i in 1: length(Dir))
  {
    T12 <- read.csv(paste0(Dir[i], "Results/", "SmCCA", "CVmeanDeltaCors.csv"))
    #pen <- which(T12[ , 4] == min(T12[ , 4]))
    ind <- which(T12[ , 5] == min(T12[ , 5]))
    temp <- T12[ind,-1]
    error <- rbind(error, as.numeric(temp))
    print(temp)
    ### Getting SD
    sdtrain_temp <- c()
    sdtest_temp <- c()
    sderror_temp <- c()
    for (j in 1:5)
    {
      error_dat <- read.csv(paste0(Dir[j], "CV_", j, "/SmCCA/SCCA_",subSamp,"_allDeltaCor.csv"))
      sdtrain_temp <- c(sdtrain_temp, error_dat[ind,4])
      sdtest_temp <- c(sdtest_temp, error_dat[ind,5])
      sderror_temp <- c(sderror_temp, error_dat[ind,6])
    }
    sd_train[i] <- sd(sdtrain_temp)
    sd_test[i] <- sd(sdtest_temp)
    sd_error[i] <- sd(sderror_temp)
  }
  error <- error[-1,]
  error <- data.frame(pen1 = error[,1], pen2 = error[,2], testcc = error[,3], pred_error = error[,4], sd_train, sd_test, sd_error)
  #error <- as.data.frame(error)
  return(error)
}

#' Grid search evaluation method that output the final scaling factor
#' @param Error The error data output from get_error() function
#' @example evaluation(my_error)
#' @export

evaluation <- function(Error)
{

  Error <- cbind(Error, name = 1:dim(Error)[1])
  index <- which.max(Error[,4])


  candidate_1 <- which((Error[index,4]- Error[1:dim(Error)[1],4])/Error[index,4]<= 0.05)

  new_error <- Error[candidate_1,]
  error_index <- which.min(new_error[,5])

  candidate_2 <- which((new_error[error_index,5]- new_error[1:dim(new_error)[1],5])/new_error[error_index,5]<= 0.05)
  decision <- new_error[which.min(new_error[candidate_2,]$name),]
  return(decision)

}


#'Enrichment analysis for network module interpretation
#'
#'This is the enrichment analysis that can be done through either fisher exact test or chi-square test
#'by changing the background gene set into the one that we specify, and obtain the enrichment pathways
#'with significance either through KEGG or GO
#'
#'@param genelist The list of gene that is used for enrichment analysis. Needs to be entrez id in character form
#'@param background The background gene set user specify, it should be in the same form as the genelist
#'@param db Database for enrichment, either KEGG or GO
#'@param species The species that experiment data is dealing with
#'@param Num At least what number of genes should be in a specific pathway
#'@param test_method Either can be Fisher's exact test or Chi-square test
#'@param adjustment P-value adjustment method, default set to FDR
#'
#'@example enrichment_GO_421_multimir <- enrichment(genelist = gene_421_multimir, test_method = "fisher", background = my_background_gene, db = "GO", species = "mmu", Num = 2, adjustment  = "fdr")
#'@export
#'
enrichment <- function(genelist,background, db, species, Num = 2, test_method = "fisher", adjustment  = "fdr")
{
  if (db == "kegg")
  {
    kegg_rno <- get_kegg(species)
    list_enrich <- diffEnrich::pathEnrich(gk_obj = kegg_rno, gene_list = genelist, cutoff = 1, N = Num)$enrich_table
    background_enrich <- diffEnrich::pathEnrich(gk_obj = kegg_rno, gene_list = background, cutoff = 1, N = Num)$enrich_table
    common_path <- as.character(intersect(list_enrich$KEGG_PATHWAY_ID, background_enrich$KEGG_PATHWAY_ID))

    ##### Creating 2 by 2 table
    Num_list <- length(genelist)
    Num_background <- length(background)
    p_value <- c()
    description <- c()
    Num_in_list <- c()
    Num_in_background <- c()
    p_value_chisq <- c()
    common_pathway <- c()
    odds_ratio <- c()

    ##### Perform the Fisher exact test or Chi-square test
    for (i in 1:length(common_path))
    {
      both <- list_enrich$KEGG_PATHWAY_in_list[list_enrich$KEGG_PATHWAY_ID == common_path[i]]
      only_background <- background_enrich$KEGG_PATHWAY_in_list[background_enrich$KEGG_PATHWAY_ID == common_path[i]]
      ct = as.table( rbind(c(both, Num_list - both), c(only_background - both, Num_background - both - Num_list - only_background)))
      or <- (ct[1,1] * ct[2,2])/(ct[1,2] * ct[2,1])
      fisher_result <- fisher.test(ct)
      chisq_result <- chisq.test(ct)
      if ((ct[1,1]/Num_list) >= (ct[2,1]/(Num_background - Num_list)))
      {
        p_value <- c(p_value, fisher_result$p.value)
        odds_ratio <- c(odds_ratio, or)
        p_value_chisq <- c(p_value_chisq, chisq_result$p.value)
        description <- c(description,list_enrich$KEGG_PATHWAY_description[list_enrich$KEGG_PATHWAY_ID == common_path[i]])
        Num_in_list <- c(Num_in_list,list_enrich$KEGG_PATHWAY_in_list[list_enrich$KEGG_PATHWAY_ID == common_path[i]])
        Num_in_background <- c(Num_in_background,background_enrich$KEGG_PATHWAY_in_list[background_enrich$KEGG_PATHWAY_ID == common_path[i]])
        common_pathway <- c(common_pathway, common_path[i])
      }
    }

  }




  if (db == "GO")
  {
    list_enrich <- clusterProfiler::enrichGO(genelist, 'org.Mm.eg.db', ont="ALL", pvalueCutoff=1, readable = TRUE)@result
    background_enrich <- clusterProfiler::enrichGO(background, 'org.Mm.eg.db', ont="ALL", pvalueCutoff=1, readable = TRUE)@result
    common_path <- as.character(intersect(list_enrich$ID, background_enrich$ID))

    ##### Creating 2 by 2 table
    Num_list <- length(genelist)
    Num_background <- length(background)
    p_value <- c()
    description <- c()
    Num_in_list <- c()
    Num_in_background <- c()
    p_value_chisq <- c()
    common_pathway <- c()
    odds_ratio <- c()

    ##### Perform the Fisher exact test or Chi-square test
    for (i in 1:length(common_path))
    {
      both <- list_enrich$Count[list_enrich$ID == common_path[i]]
      only_background <- background_enrich$Count[background_enrich$ID == common_path[i]]
      ct = as.table( rbind(c(both, Num_list - both), c(only_background - both, Num_background - both - Num_list - only_background)))
      or<- (ct[1,1] * ct[2,2])/(ct[1,2] * ct[2,1])
      if (dim(ct)[1] == 0|dim(ct)[2] == 0)
      {stop("This enrichment analysis is not returning any result")}

      fisher_result <- fisher.test(ct)
      chisq_result <- chisq.test(ct)
      if ((ct[1,1]/Num_list) >= (ct[2,1]/(Num_background - Num_list)))
      {
        odds_ratio <- c(odds_ratio, or)
        p_value <- c(p_value, fisher_result$p.value)
        p_value_chisq <- c(p_value_chisq, chisq_result$p.value)
        description <- c(description,list_enrich$Description[list_enrich$ID == common_path[i]])
        Num_in_list <- c(Num_in_list,list_enrich$Count[list_enrich$ID == common_path[i]])
        Num_in_background <- c(Num_in_background,background_enrich$Count[background_enrich$ID == common_path[i]])
        common_pathway <- c(common_pathway, common_path[i])
      }
    }



  }


  if (test_method == "fisher")
  {
    p_adjust <- p.adjust(p_value, method = adjustment)
    output_result <- data.frame(pathway = common_pathway, description = description, odds_ratio = odds_ratio, p_value = p_value, p_adjust = p_adjust, Num_in_list = Num_in_list, Num_in_background = Num_in_background)
    output_result <-output_result[order(output_result$p_adjust),]
  }
  if (test_method == "chisq") {
    p_adjust <- p.adjust(p_value_chisq, method = adjustment)
    output_result <- data.frame(pathway = common_pathway, description = description, odds_ratio = odds_ratio, p_value = p_value_chisq, p_adjust = p_adjust, Num_in_list = Num_in_list, Num_in_background = Num_in_background)
    output_result <-output_result[order(output_result$p_adjust),]
  }
  if ((test_method != "fisher") & (test_method != "chisq")) {stop('The test method being specified is wrong, try fisher or chisq.')}

  cat("----------------------------------------------------\n")
  cat("Output the Top Pathways From the Enrichment Analysis\n")
  cat("----------------------------------------------------\n")
  print(output_result[1:10,])
  if (output_result$p_adjust[1] > 0.05){
    warning("This enrichment analysis did not provide significant result")}
  print(output_result)
}


#'Enrichment analysis for genes in the network modules
#'
#'The function enrichment() is dealing with the multimiR target genes, enrich network here is dealing with
#'network genes for each network module
#'@param gene_name Gene symbol for the module, in a form of character vector
#'@param background Background gene set user specified
#'@param db Database of interest: KEGG or GO
#'
#'@example gene_421 <- read.csv("C:/Research Material/formalresult/enrichment_table/gene_folder/GeneID_Net_421.csv", header = FALSE)
#'my_421 <- unlist(mapIds(org.Mm.eg.db, gene_421_multimir,  'SYMBOL','ENTREZID'))
#'intersect(my_421, as.character(gene_421$V1))
#'enrich_network_421 <- Enrich_Networks(as.character(gene_421$V1), gene_background, db = "kegg")
#'@export

Enrich_Networks <- function(gene_name, background, db)
{
  genelist <- select(org.Mm.eg.db, gene_name, "ENTREZID","SYMBOL")
  genelist <- genelist[!is.na(genelist$ENTREZID),]
  genelist <- as.character(genelist$ENTREZID)
  if (db == "kegg")
  {
    enrichment_result <- enrichment(genelist, test_method = "fisher", background, db = "kegg", species = "mmu", Num = 1, adjustment  = "fdr")
  }
  if (db == "GO")
  {

    enrichment_result <- enrichment(genelist, test_method = "fisher", background = as.character(background), db = "GO", species = "mmu", Num = 1, adjustment  = "fdr")
  }
  return(enrichment_result)
}








