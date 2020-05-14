#' Running SmCCNet with criterion specified by users
#'
#' @param dir the working directory
#' @param X1 first dataset
#' @param X2 second dataset
#' @param Y  phenotype
#' @param factor_1  scaling factor 1
#' @param factor_2  scaling factor 2
#' @param factor_3  scaling factor 3
#' @param pheno     name of phenotype
#' @param K         K-folds cross-validation
#' @param s1        subsampling parameter 1
#' @param s2        subsampling parameter 2
#' @param subSamp   subsampling times
#' @param pen1      first penalty sequence
#' @param pen2      second penalty sequence
#' @param Type      parallel processing type
#' @return SmCCNet result stored in the local directory
#' @example library(SmCCNet)
#' X1 <- geneExpr[,1:100]
#' X2 <- mirnaExpr[,1:50]
#' Y <- pheno
#' dim(X1)
#' dim(X2)
#' relative <- 1:20
#' factor_1 <- 3/(1+relative *2)
#' factor_2 <- factor_3 <- (3 - factor_1)/2
#' CVDirs <- SmCCNet_Running(dir = "C:/Research Material/formalresult/",X1,X2,Y,factor_1,factor_2,factor_3, pheno = "Pheno", K = 5,s1 = 0.7, s2 = 0.9, subSamp = 10,
#' pen1 = seq(.05, .3, by = .025),
#' pen2 = seq(.05, .3, by = .025), Type = "PSOCK")
#' @importFrom SmCCNet getRobustPseudoWeights
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterEvalQ
#' @export

SmCCNet_Running <- function(dir,X1,X2,Y,factor_1,factor_2,factor_3, pheno, K,s1 = 0.7, s2 = 0.9, subSamp = 1000,
                            pen1 = seq(.05, .3, by = .025),
                            pen2 = seq(.05, .3, by = .025), Type = "PSOCK")
{
  #utils::globalVariables(c("AbarLabel", "p", "CVDirs", "CVDir"))
  dir.create(paste0(dir,pheno))
  X <- cbind(X1, X2)
  bigCor2 <- cor(X)
  AbarLabel <- colnames(bigCor2)
  p <- ncol(bigCor2)
  CVDirs <- rep(NA,20)
  for (i in 1:length(factor_1))
  {
    CVDirs[i] <- paste0(dir, pheno,  "/", K, "foldCV", pheno, "_CCcoef_1_",
                         factor_1[i], "_", factor_2[i], "/")
    dir.create(CVDirs[i])
  }
  setwd(dir)
  ########### Running SmCCNet
  for (i in 1:length(factor_1))
  {
    #library(parallel)
    #library(SmCCNet)
    CCcoef <- c(factor_1[i], factor_2[i], factor_3[i])
    CVDir <- CVDirs[i]
    set.seed(12345)

    N <- nrow(X1); p1 <- ncol(X1); p2 <- ncol(X2)
    #s1 <- 0.7; s2 <- 0.9; subSamp <- 1000
    #pen1 <- seq(.05, .3, by = .025)
    #pen2 <- seq(.05, .3, by = .025)
    P1P2 <- expand.grid(pen1, pen2)
    a <- sqrt(p1 * s1) * P1P2[ , 1]; a[a < 1] <- 1
    b <- sqrt(p2 * s2) * P1P2[ , 2]; b[b < 1] <- 1
    P1P2 <- P1P2[which(a>b), ]
    save(X1, X2, Y, s1, s2, subSamp, P1P2, CCcoef,
         file = paste0(CVDir, "Data.Rdata"))

    foldIdx <- split(1:N, sample(1:N, K))
    for(i in 1:K){
      iIdx <- foldIdx[[i]]
      x1.train <- scale(X1[-iIdx, ])
      x2.train <- scale(X2[-iIdx, ])
      yy.train <- scale(Y[-iIdx, ])
      x1.test <- scale(X1[iIdx, ])
      x2.test <- scale(X2[iIdx, ])
      yy.test <- scale(Y[iIdx, ])

      if(is.na(min(min(x1.train), min(x2.train), min(yy.train), min(x1.test), min(x2.test), min(yy.test)))){
        stop("Invalid scaled data.")
      }

      subD <- paste0(CVDir, "CV_", i, "/")
      dir.create(subD)
      save(x1.train, x2.train, yy.train, x1.test, x2.test, yy.test,
           s1, s2, P1P2, p1, p2, subSamp, CCcoef,
           file = paste0(subD, "Data.Rdata"))
    }



    # # Note that there are 8 cores on my machine.
    # # detectCores()
    # # [1] 8
    #
    # # Initiate cluster. Use type "FORK" to automatically include all environment variables.
    # cl <- makeCluster(detectCores() - 1, type = "FORK")
    cl <- makeCluster(5, type = Type)
    clusterEvalQ(cl, library(SmCCNet))
    #clusterEvalQ(cl, CVDir <- CVDirs[i])
    clusterExport(cl = cl, "CVDir", envir = environment())
    parSapply(cl, 1:5, function(CVidx){
      # source("../../../ModifiedPMA.R")
      # source("../../../SmCCNetSource.R")

      subD <- paste0(CVDir, "CV_", CVidx, "/")
      load(paste0(subD, "Data.Rdata"))
      dir.create(paste0(subD, "SmCCA/"))

      RhoTrain <- RhoTest <- DeltaCor <- rep(0, nrow(P1P2))
      for(idx in 1:nrow(P1P2)){
        # print(paste0("Running SmCCA on CV_", CVidx, " idx=", idx))

        l1 <- P1P2[idx, 1]; l2 <- P1P2[idx, 2]
        Ws <- getRobustPseudoWeights(x1.train, x2.train, yy.train, l1, l2, s1, s2,
                                     NoTrait = FALSE, FilterByTrait = FALSE,
                                     SubsamplingNum = subSamp, CCcoef = CCcoef)
        # Ws <- Matrix(Ws)
        meanW <- rowMeans(Ws)
        v <- meanW[1:p1]; u <- meanW[p1 + 1:p2]

        rho.train <- cor(x1.train %*% v, x2.train %*% u) * CCcoef[1] +
          cor(x1.train %*% v, yy.train) * CCcoef[2] +
          cor(x2.train %*% u, yy.train) * CCcoef[3]
        rho.test <- cor(x1.test %*% v, x2.test %*% u) * CCcoef[1] +
          cor(x1.test %*% v, yy.test) * CCcoef[2] +
          cor(x2.test %*% u, yy.test) * CCcoef[3]

        RhoTrain[idx] <- round(rho.train, digits = 5)
        RhoTest[idx] <- round(rho.test, digits = 5)
        DeltaCor[idx] <- abs(rho.train - rho.test)

        if(idx %% 10 == 0){
          save(P1P2, RhoTrain, RhoTest, DeltaCor, idx,
               file = paste0(subD, "temp.Rdata"))
        }

      }

      DeltaCor.all <- cbind(P1P2, RhoTrain, RhoTest, DeltaCor)
      colnames(DeltaCor.all) <- c("l1", "l2", "Training CC", "Test CC", "CC Pred. Error")
      write.csv(DeltaCor.all,
                file = paste0(subD, "SmCCA/SCCA_", subSamp,"_allDeltaCor.csv"))

      system(paste0("rm ", subD, "temp.Rdata"))

      return(CVidx)

    }
    )



    # Close cluster
    stopCluster(cl)
  }
  return(CVDirs)
}


#' Running overall SmCCNet based on previous cross-validation result
#' @param dir working directory
#' @param CVDirs result directory
#' @param X1 first dataset
#' @param X2 second dataset
#' @param Y phenotype
#' @param factor_1 first scaling factor
#' @param factor_2 second scaling factor
#' @param factor_3 third scaling factor
#' @param pheno phenotype name
#' @param cutting a sequence of cutting edge
#' @param subSamp number of subsampling
#' @param s1 first subsampling parameter
#' @param s2 second subsampling parameter
#' @param K cross-validation parameter
#' @param phenos name of phenotype
#' @example overall_running(dir = "C:/Research Material/formalresult/",CVDirs, X1,X2,Y,factor_1,factor_2,factor_3, pheno = "Pheno", cutting = c(0,0.01,0.05, 0.3),subSamp = 10,
#' s1 = 0.7, s2 = 0.9, K = 5, phenos = "Pheno")
#' @return SmCCNet result stored in the local directory
#' @export

overall_running <-function(dir,CVDirs, X1,X2,Y,factor_1,factor_2,factor_3, pheno, cutting = c(0,0.01,0.3), subSamp = 1000,
                           s1 = 0.7, s2 = 0.9, K = 5, phenos)
{
  setwd(dir)
  #CVDirs <- rep(NA,20)
  #for (i in 1:length(factor_1))
  #{
    #CVDirs[i] <- paste0(dir,"/", pheno,  "/", K, "foldCV", pheno, "_CCcoef_1_",
                        #factor_1[i], "_", factor_2[i], "/")
    #dir.create(CVDirs[i])
  #}
  source("SCCAdiagTools.R")
  X <- cbind(X1, X2)
  bigCor2 <- cor(X)
  AbarLabel <- colnames(bigCor2)
  p <- ncol(bigCor2)
  N <- nrow(X1); p1 <- ncol(X1); p2 <- ncol(X2)
  for(i in 1:length(factor_1)){
    CCcoef <- c(factor_1[i], factor_2[i], factor_3[i])
    CVDir <- CVDirs[i]
    plotD <- paste0(CVDir, "Figures/")
    saveD <- paste0(CVDir, "Results/")
    # # Create contour plots for canonical correlation prediction error.
    plotCVcontour(CVDir, "SmCCA", NumSubsamp = subSamp)


    # Run SCCA on the entired dataset based on CV proposed penalty parameters.
    for(Method in "SmCCA"){
      T12 <- read.csv(paste0(CVDir, "Results/", Method, "CVmeanDeltaCors.csv"))[ , -1]
      pen <- which(T12[ , 4] == min(T12[ , 4]))

      if(Method == "SmCCA"){
        FilterByTrait <- FALSE
      }else if(Method == "SsCCA"){
        FilterByTrait <- TRUE
      }

      l1 <- T12$l1[pen]; l2 <- T12$l2[pen]
      Ws <- getRobustPseudoWeights(X1, X2, Y, l1, l2, s1, s2, NoTrait = FALSE,
                                   FilterByTrait = FilterByTrait,
                                   SubsamplingNum = subSamp, CCcoef = CCcoef)
      save(l1, l2, X1, X2, Y, s1, s2, Ws, CCcoef,
           file = paste0(saveD, Method, K, "foldSamp", subSamp, "_", pen, ".Rdata"))

      Abar <- getAbar(Ws, FeatureLabel = AbarLabel[1:p])
      save(l1, l2, X1, X2, Y, s1, s2, Ws, Abar, CCcoef,
           file = paste0(saveD, Method, K, "foldSamp", subSamp, "_", pen, ".Rdata"))


      sccafile <- paste0(saveD, Method, K, "foldSamp", subSamp, "_", pen, ".Rdata")
      saveresult <- paste0(saveD, Method, K, "foldSamp", subSamp, "_", pen, "Result.Rdata")
      penDir <- paste0(plotD, Method, "pen", pen, "/")
      dir.create(penDir)

      mirGeneModule <- getMultiOmicsModules(Abar, p1, PlotTree = FALSE)
      save(pen, penDir, l1, l2, mirGeneModule, file = saveresult)

      for(edgeCut in cutting){
        if(edgeCut == 0){
          netD <- paste0(penDir, "NetworkEdgeCutpt0/")
          titlePre <- paste0("Net ")
          savePre <- paste0(netD, "Net")
        }else{
          netD <- paste0(penDir, "NetworkEdgeCutpt", substring(edgeCut, 3), "/")
          titlePre <- paste0("Trimmed Net ")
          savePre <- paste0(netD, "TrimmedNet")
        }
        dir.create(netD)

        lapply(1:length(mirGeneModule), function(z){
          saveplot.z <- paste0(savePre, z, ".pdf")
          plottitle <- paste0(titlePre, z)
          plotMultiOmicsNetwork(Abar, bigCor2, mirGeneModule, z, P1 = p1,
                                EdgeCut = edgeCut, FeatureLabel = AbarLabel,
                                AddCorrSign = TRUE, SaveFile = saveplot.z,
                                ShowType1Label = TRUE, ShowType2Label = TRUE,
                                PlotTitle = plottitle, NetLayout = "circle",
                                ShowNodes = TRUE, VertexLabelCex = 1, VertexSize = 1)
        })
      }
    }
  }
}

#' Extracting Omics modules, heatmap and other grid search related measures
#' @param CVDirs result directory based on the previous SmCCNet running directory, which
#' is given by the output of overall_running() function
#' @param grids this is the number of grids that are used for the grid search, which is given by length(factor_1)
#' @param phenos the name of the phenotype
#' @param cutting a vector of cutting edge we are using, which should be consistent with overall_running()
#' @example Extract_Information(CVDirs, grid = length(factor_1), phenos = 'Pheno', cutting = c(0,0.01,0.05,0.3))
#' @return summary results in the result directory
#' @export


Extract_Information <- function(CVDirs, grids, phenos, cutting)
{

  X <- cbind(X1, X2)
  bigCor2 <- cor(X)
  AbarLabel <- colnames(bigCor2)
  p <- ncol(bigCor2)
  N <- nrow(X1); p1 <<- ncol(X1); p2 <<- ncol(X2)
Pens <- Abars <- AbarSigneds <- mirGeneModules <- vector("list", 20)
for(ind in (1:length(factor_1))){
    #pheno <- phenos[2]
    sccaFs <- paste0(CVDirs[ind], "Results/SmCCA*.Rdata")
    sapply(Sys.glob(sccaFs), load, .GlobalEnv)
    Pens[[ind]] <- pen; mirGeneModules[[ind]] <- mirGeneModule
    Abars[[ind]] <- Abar
    colnames(Abars[[ind]]) <- rownames(Abars[[ind]]) <- AbarLabel
    AbarSigneds[[ind]] <- Abar * sign(bigCor2)
  }

  for(j in 1:length(cutting))
  {
  corr_pca <- list()
  varexp_pca <- list()
  final_result_pca <- list()
  variance_exp_pca <-c()
  modulemirna_pca <- c()
  modulemrna_pca <- c()

  for (i in 1:length(factor_1))
  {
    mirGeneM <- mirGeneModules[[i]]
    Abar <- Abars[[i]]
    cvDir <- CVDirs[i]
    dir.create(paste0(cvDir,"indvCor"))
    dir.create(paste0(cvDir,"ModulePC1"))
    cutMirGeneM <- reducedMirGeneModule(mirGeneM, Abar, edgeCut = cutting[j])
    cutpt <- paste0("CutPt",cutting[j] * 100)
    result <- checkCorr(cutMirGeneM, Abar, X, Y, phenos, cvDir, cutpt, ycoord = NULL)
    corr_pca[[i]] <- result[[1]]
    varexp_pca[[i]] <- result[[2]]
    variance_exp_pca <- c(variance_exp_pca, result[[3]])
    modulemirna_pca <- c(modulemirna_pca, result[[4]])
    modulemrna_pca <- c(modulemrna_pca, result[[5]])
    names(corr_pca)[[i]] <- paste0("1:",i,":",i, "correlation")
    names(varexp_pca)[[i]] <- paste0("1:",i,":",i, "variance")
  }
  final_result_pca <- c(corr_pca, varexp_pca)

  mirGeneModules.trimmed <- sapply(1:20, function(x){
    y <- reducedMirGeneModule(mirGeneModules[[x]], Abars[[x]], edgeCut = cutting[j])
    y <- Filter(Negate(is.null), y)
    return(y)
  })

  ############################################### Gene Extraction

  geneLabel <- AbarLabel[1:p1]
  mirLabel <-  AbarLabel[p1 + 1:p2]
  geneLists <- mirLists <- vector("list", 20)
  for(ind in 1:length(factor_1)){
    #pheno <- phenos[ind]
    for (I in 1:length(mirGeneModules.trimmed[[ind]])){

    cvDir <- CVDirs[ind]
    mirGeneM <- mirGeneModules.trimmed[[ind]][I]

    #dir.create(cvDir, "MirGeneID")
    dir.create(paste0(cvDir, "MirGeneID"))
    saveD <<- paste0(cvDir, "MirGeneID/CutPt",cutting[j],"_",I,"/")
    dir.create(saveD)

    geneLists[[ind]] <- sapply(1:length(mirGeneM), function(x){
      saveF <- paste0(saveD, "GeneID_Net", ".csv")
      genes <- getGeneID(mirGeneM, x, geneLabel, saveF)
      return(genes)
    })
    mirLists[[ind]] <- sapply(1:length(mirGeneM), function(x){
      saveF <- paste0(saveD, "MirID_Net", ".csv")
      genes <- getMirID(mirGeneM, x, p1, mirLabel, saveF)
      return(genes)
    })
    }
  }

}
for (i in 1:length(factor_1))
{
  extract_heatmap(CVDirs[i], cutting)
}

}




############################### Source Function
#' (INTERNAL)
#' @keywords internal


reducedMirGeneModule <- function(mirGeneM, abar, edgeCut = edgeCut){
  # Reduce mirGeneModule list according to the edge cut provided.

  reducedMirGeneM <- lapply(mirGeneM, function(x){
    subAbar <- abar[x, x]
    edgeMax <- apply(subAbar, 1, max)
    m <- which(edgeMax >= edgeCut)
    if(length(m) > 0){y <- x[m]}else(y <- NULL)
    return(y)
  })

  return(reducedMirGeneM)
}




#' (INTERNAL)
#' @keywords internal

checkCorr <- function(cutModule, abar, X, Y, phenos, cvdir, cutname, ycoord = NULL){
  # Find module trait relationships through eigengenes
  moduleColors <- rep(0, nrow(abar))
  K1 <- length(cutModule)

  for(k in 1:K1){
    M <- cutModule[[k]]
    if(!is.null(M)){moduleColors[M] <- k}
  }
  MEs <- moduleEigengenes(X, moduleColors)$eigengenes
  moduleTraitCor <- abs(cor(MEs, Y))
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(X))
  propVar <- propVarExplained(X, moduleColors, MEs, corFnc = "cor")
  textMatrix <- paste0(signif(moduleTraitCor, 2), " (",
                       signif(moduleTraitPvalue, 1), ")");
  dim(textMatrix) <- dim(moduleTraitCor)
  heatmap_info <- data.frame(moduleTraitCor, colnames(MEs), textMatrix)
  write.csv(heatmap_info, file = paste0(cvdir, "ModulePC1/heatmap_info",
                                        cutname, ".csv"))

  # Create a table that saves module information
  cutModule2 <- Filter(Negate(is.null), cutModule)
  K2 <- length(cutModule2)
  mirNum <- geneNum <- rep(0, K2)
  for(k in 1:K2){
    M <- cutModule2[[k]]
    mirNum[k] <- sum(M > p1)
    geneNum[k] <- sum(M <= p1)
  }
  mirNum <- c(p2 - sum(mirNum), mirNum)
  geneNum <- c(p1 - sum(geneNum), geneNum)
  moduleT <- cbind(signif(moduleTraitCor, 4), signif(moduleTraitPvalue, 2),
                   signif(propVar, 3),
                   mirNum, geneNum, mirNum+geneNum)
  rownames(moduleT) <- colnames(MEs)
  colnames(moduleT) <- c("Corr", "p-val"   , "propVarExp"
                         , "mirNum",
                         "geneNum", "modSize")
  write.csv(moduleT, file = paste0(cvdir, "ModulePC1/ModulePhenoCorpca",
                                   cutname, ".csv"))


  # Plot univariate correlation to pheno for mRNA-miRNA module members
  type <- c(rep("Gene", p1), rep("miRNA", p2))
  CorToPheno <- cor(X, Y)
  corToPheno.frame <- data.frame(CorToPheno = CorToPheno, Module = moduleColors,
                                 Type = type)

  if(cutname == "noCut" | cutname == "noCut_coarseGrid"){yLab <- "Module"}else{yLab <- "Trimmed"}
  ggp <- ggplot(corToPheno.frame, aes(x = CorToPheno, y = Module)) +
    geom_point(size = 3, alpha = 0.7, aes(shape = Type), colour = "blue") +
    ylab(yLab) + xlab ("Correlation") + xlim(c(-1,1)) +
    scale_y_reverse(limits = c(K1 , 0), breaks = 0:K1) +
    theme(axis.text = element_text(size = rel(3)),
          axis.title = element_text(size = rel(3), face = "bold")) +
    theme(legend.position = "none")

  if(is.null(ycoord)){
    ggp
  }else{
    ggp + coord_cartesian(ylim = c(0, ycoord))
  }

  ggsave(paste0(cvdir, "indvCor/IndvCorTopca", phenos, cutname, ".pdf"),
         width = 10, height =2.7)
  return(list(moduleTraitCor, propVar, propVar, mirNum, geneNum))
}



#' (INTERNAL)
#' @keywords internal

extract_heatmap <- function (data_source, cut)
{
  heatmap_data <- data.frame()
  for (i in 1:length(cut))
  {
    dat <- read.csv(paste0(data_source, "ModulePC1/heatmap_infoCutPt",
                           cut[i]*100, ".csv"))
    dat <- dat[,-1]
    ### round correlation
    dat[,1] <- round(dat[,1], digits = 2)
    if(i == 1)
    {
      mod <- 0:max(as.numeric(gsub("ME", "", dat$colnames.MEs. )))
      #PC_null <- rep(0, max(mod) + 1)
      modules_num <- dat$colnames.MEs.
      heatmap_data <- data.frame(mod)
    }
    dat$mod <- as.numeric(gsub("ME", "", dat$colnames.MEs. ))
    dat$colnames.MEs. <- NULL
    colnames(dat)[c(1,2)] <- c(paste0("PC_cut", i), paste0("text_cut", i))
    heatmap_data <- merge(heatmap_data,dat, by.x = "mod", all.x = TRUE)
    #is.na(heatmap_data[,2*i]) <- 0
    #is.na(heatmap_data[,1+2*i]) <- "NA"
  }
  heatmap_data[is.na(heatmap_data)] <- round(0,digits = 1)
  module_correlation <- heatmap_data[,grepl( "PC" , names( heatmap_data ) )]
  text_mat <- heatmap_data[,grepl( "text" , names( heatmap_data ) )]
  pdf(paste0(data_source,'indvCor/',"heatmap.pdf"),width = 20, height = 10)
  labeledHeatmap(Matrix = as.matrix(module_correlation), xLabelsAngle = 0, xLabelsAdj = 0.5,
                 xLabels = cut, cex.lab.x = 2.5,cex.lab.y = 2.5,
                 yLabels = modules_num, ySymbols = modules_num, cex.main = 2.5,
                 colorLabels = FALSE, colors = blueWhiteRed(50),
                 textMatrix = as.matrix(text_mat), textAdj = c(0.25, 0.25),#w/ corr and pvalue in cells
                 cex.text = 2.5, setStdMargins = T, zlim = c(-1,1),
                 main = "", plotLegend = FALSE
  )
  dev.off()
}





#' (INTERNAL)
#' @keywords internal

getGeneID <- function(MirGeneModule, Net, GeneLabel, SaveFile){
  # Extract gene names from a subnetwork.
  #
  # MirGeneModule: mirGeneModule list.
  # Net: number of subnetwork.
  # GeneLabel: gene names.
  # SafeFile: file name for the extract gene names.

  p1 <- length(GeneLabel)
  module <- MirGeneModule[[Net]]
  genes <- sort(as.character(GeneLabel[module[which(module <= p1)]]))
  write.table(genes, quote = FALSE, file = SaveFile, row.names = FALSE,
              col.names = FALSE)
  return(genes)
}

#' (INTERNAL)
#' @keywords internal

getMirID <- function(MirGeneModule, Net, P1, MirLabel, SaveFile){
  # Extract gene names from a subnetwork.
  #
  # MirGeneModule: mirGeneModule list.
  # Net: number of subnetwork.
  # P1: number of genes
  # MirLabel: miRNA names.
  # SafeFile: file name for the extract gene names.

  module <- MirGeneModule[[Net]]
  mirs <- sort(as.character(MirLabel[module[which(module > p1)] - p1]))
  write.table(mirs, quote = FALSE, file = SaveFile, row.names = FALSE,
              col.names = FALSE)
  return(mirs)
}



