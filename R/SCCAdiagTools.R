##################################
# Tools for getting SCCA results #
##################################
library(WGCNA)
library(igraph)
library(plotly)
library(webshot)
library(data.table)
allowWGCNAThreads()

#' (INTERNAL)
#' @keywords internal

PlotEntropyContour <- function(HTable, SavePlotFile = NULL){
  hmelt <- melt(HTable, id.vars = c("l1", "l2"))

  p <- plot_ly(hmelt, x = ~l1, y = ~l2, z = ~value, type = "contour")
  if(!is.null(SavePlotFile)){
    export(p, file = SavePlotFile)
  }
  p
}

### Obtain cross validation (CV) results.
#' (INTERNAL)
#' @keywords internal

plotCVcontour <- function(CVDir, SCCAmethod, K = 5, NumSubsamp = 1000){
    plotD <- paste0(CVDir, "Figures/")
    saveD <- paste0(CVDir, "Results/")
    dir.create(plotD); dir.create(saveD)

    testCC <- predError <- NULL
    for(j in 1:K){
        resultT <- paste0(CVDir, "CV_", j, "/", SCCAmethod,
        "/SCCA_", NumSubsamp, "_allDeltaCor.csv")
        dCorT <- read.csv(resultT)[ , -1]
        testCC <- cbind(testCC, abs(dCorT[ , 4]))
        predError <- cbind(predError, dCorT[ , 5])
    }

    S1 <- rowMeans(testCC)
    S2 <- rowMeans(predError)

    T12 <- dCorT[ , -3]; T12[ , 3] <- S1; T12[ , 4] <- S2
    write.csv(T12, file = paste0(saveD, SCCAmethod, "CVmeanDeltaCors.csv"))

    print(paste0("testCC choice: ", which(S1 == min(S1))))
    print(paste0("CC Pred. Error choice: ", which(S2 == min(S2))))

    PlotEntropyContour(T12[ , -4], paste0(plotD, SCCAmethod, "meanTestCCContour.pdf"))
    PlotEntropyContour(T12[ , -3], paste0(plotD, SCCAmethod, "meanDeltaCorContour.pdf"))
}

#' (INTERNAL)
#' @keywords internal

GetBestIdx <- function(DirName, csvFile = TRUE){
  # Report the index for the best SCCA result.
  #
  # DirName: directory where the data was saved in.
  # csvFile: a table file named "SCCA_100_topH.csv" that saves entropy values
  #   for each result file.

  if(csvFile){
    H.tab <- read.csv(paste0(DirName, "SCCA_100_topH.csv"))
    idx <-  H.tab$X[which(H.tab$H == min(H.tab$H))]
  }else{
    sccaFiles <- Sys.glob(paste0(DirName, "Bipartite/*"))
    Hs <- sapply(sccaFiles, function(x){load(x); return(H)})
    nchar.dir <- nchar(DirName)
    idces <- sapply(sccaFiles, function(x){
      nchar.x <- nchar(x)
      idx.x <- substr(x, nchar.dir+20, nchar.x-6)
      return(idx.x)
    })
    idx <- idces[which(Hs == min(Hs))]
  }
  return(idx)
}

#' (INTERNAL)
#' @keywords internal

PlotAbar <- function(sccaFile, SavePlotFile, Zoom = NULL){
  # Plot the average weight matrix.
  #
  # sccaFile: SCCA result file.
  # SavePlotFile: file name for the plot to be saved as.
  # Zoom: zoom-in vector.

  load(sccaFile)
  pdf(SavePlotFile)
  if(is.null(Zoom)){zoom <- 1:nrow(Abar)}else{zoom = Zoom}
  A <- matrix(t(Abar[zoom, zoom]), nrow = length(zoom))
  image(zoom, zoom, A, col = grey(seq(1, 0, length = 2560)),
        xlab = "", ylab = "",
        main = paste0("Edge weight matrix (pen=(", l1, ",", l2,"))"))
  dev.off()
}

#' (INTERNAL)
#' @keywords internal

GetFuncGrp <- function(sccaFile, SaveResultFile, SavePlotFile = NULL,
                       Cutoff = 1,
                       GetPerformance = FALSE, Signal = NULL){
  # Report functional group obtained using tree cut.
  #
  # sccaFile: SCCA result file.
  # SaveResultFile: file name for the result to be saved as.
  # SavePlotFile: file name for the hierarchical tree plot to be saved as.
  # Cutoff: tree cut height.
  # GetPerformance: whether to compute precision and recall.
  # Signal: indeces for the signals.


  load(sccaFile)
  rownames(Abar) <- colnames(Abar) <- 1:nrow(Abar)
  hc <- hclust(as.dist(1 - Abar))

  if(is.null(SavePlotFile)){
    plot(hc, main = paste0("pen=(", l1, ",", l2,")"), cex = .6)
  }else{
    pdf(width = 16, SavePlotFile)
    plot(hc, main = paste0("pen=(", l1, ",", l2,")"), cex = .6)
    dev.off()
  }

  cut.merge <- hc$merge[hc$height < Cutoff, ]
  lower.leaves <- sort(-cut.merge[cut.merge<0])
  precision <- recall <- NA
  if(GetPerformance){
    TP <- length(intersect(lower.leaves, Signal))
    FP <- length(setdiff(lower.leaves, Signal))
    FN <- length(setdiff(Signal, lower.leaves))
    precision <- TP/(TP + FP)
    recall <- TP/(TP + FN)
    print(list(precision = precision, recall = recall, signal = Signal))
  }
  save(hc, l1, l2, Cutoff, lower.leaves, Signal, precision, recall,
       file = SaveResultFile)

  return(lower.leaves)
}

#' (INTERNAL)
#' @keywords internal

Leaves2Module <- function(Leaves, GrpID){
  # Create module index list.
  #
  # Leaves: indices for the nodes below the threshold in the hierarchical tree.
  # GrpID: group membership for all nodes.

  id <- GrpID[Leaves]
  M <- lapply(1:length(unique(id)), function(x){
    M.x <- Leaves[which(id == unique(id)[x])]
    return(M.x)
  })

  return(M)
}

#' (INTERNAL)
#' @keywords internal

MakePerformanceTable <- function(ResultFiles, RowNames = NULL,
                                 SaveTableFile = NULL){
  # Consolidate performance report from several result files.
  #
  # ColNames: desired colnames for the report table.
  # ResultFiles: file names of the results to be reported.

  N <- length(ResultFiles)
  Perform <- matrix(NA, ncol = 4, nrow = N)
  for(n in 1:N){
    load(ResultFiles[n])
    Perform[n, ] <- c(l1, l2, precision, recall)
  }
  Perform <- round(Perform, digits = 4)
  colnames(Perform) <- c("l1", "l2", "Precision", "Recall")
  rownames(Perform) <- RowNames

  if(!is.null(SaveTableFile)){write.csv(Perform, SaveTableFile)}

  return(Perform)
}

#' (INTERNAL)
#' @keywords internal

PickWGCNApower <- function(X, ResultDir, PlotDir, Network = "unsigned",
                           Corr = "cor"){
  # Choose the network soft threshold.
  #
  # X: n by p data matrix.
  # ResultDir: directory where the results to be saved in.
  # PlotDir: directory where the plots to be saved in.
  # Network: type of network, either "unsigned" or "signed".
  # Corr: correlation used in adjacency calculation, either "cor" or "bicor".

  wgcnaPlot <- paste0(PlotDir, "WGCNA", Network, "/")
  wgcnaResult <- paste0(ResultDir, "WGCNA", Network, "/")
  dir.create(wgcnaPlot); dir.create(wgcnaResult)

  powers <- 1:25
  # Choose a set of soft-thresholding powers
  sft <- pickSoftThreshold(X, powerVector = powers, verbose = 3,
                          networkType = Network, corFnc = Corr)
  save(powers, sft, file = paste0(wgcnaResult, "sft0.Rdata"))
  collectGarbage();

  # Plot the results:
  # Will plot these columns of the returned scale free analysis tables
  plotCols = c(2,5,6,7)
  colNames = c("Scale Free Topology Model Fit", "Mean connectivity",
               "Median connectivity", "Max connectivity")

  # Plot the quantities in the chosen columns vs. the soft thresholding power
  pdf(file = paste0(wgcnaPlot, "sft0.pdf"), width = 12, height = 12)
  par(mfcol = c(2,2));
  par(mar = c(4.2, 4.2 , 2.2, 0.5))
  ylim = matrix(NA, nrow = 2, ncol = 4)
  for (col in 1:length(plotCols))
  {
    ylim[1, col] = min(ylim[1, col], sft$fitIndices[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], sft$fitIndices[, plotCols[col]], na.rm = TRUE);
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)", ylab=colNames[col],type="n",
         main = colNames[col], ylim = ylim[, col]);
    addGrid();
    if (col==1){
      # Scale-free topology fit index as a function of the soft-thresholding power
      text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
           labels=powers, col="blue");
      # this line corresponds to using an R^2 cut-off of h
      abline(h=0.85,col = "red", lty = 2)
    }else{
      text(sft$fitIndices[,1], sft$fitIndices[,plotCols[col]],
           labels=powers, col="blue");
    }
  }
  dev.off()

  return(sft$powerEstimate)
}

#' (INTERNAL)
#' @keywords internal

RunWGCNA <- function(X, wgcnaResult, wgcnaPlot, Power = 6,
                     Network = "unsigned", Corr = "cor",
                     minModuleSize = 3, mergeCutHeight = 0.25){
  # Run WGCNA analysis.
  #
  # X: n by p data matrix.
  # Power: network soft threshold.
  # wgcnaResult: directory where the results to be saved in.
  # wgcnaPlot: directory where the plots to be saved in.
  # Network: type of network, either "unsigned" or "signed".
  # Corr: correlation used in adjacency calculation, either "cor" or "bicor".
  # minModuleSize: minimum module size for module detection.
  # mergeCutHeight: dendrogram cut height for module merging.

  dir.create(paste0(wgcnaResult, "Beta", Power))
  if(wgcnaResult != wgcnaPlot){dir.create(paste0(wgcnaPlot, "Beta", Power))}


  #Get modules and plot the network
  net <- blockwiseModules(X, power = Power, TOMType = Network,
                          minModuleSize = minModuleSize,
                          mergeCutHeight = mergeCutHeight)

  if(Network == "unsigned"){signed = FALSE}else{signed = TRUE}
  pdf(file = paste0(wgcnaPlot, "Beta", Power, "/MENetwork.pdf"),
      width = 12, height = 12)
  plotEigengeneNetworks(net$MEs, paste0("LXS (beta=", Power, ")"),
                        signed = signed, plotDendrograms = FALSE)
  dev.off()
  pdf(file = paste0(wgcnaPlot, "Beta", Power, "/Block1MENetwork.pdf"),
      width = 12, height = 12)
  plotDendroAndColors(
    net$dendrograms[[1]], net$colors[net$blockGenes[[1]]], "Module colors",
    main = "Gene dendrogram and module colors",
    dendroLabels = NULL, addGuide = TRUE)
  dev.off()

  save(X, net, file = paste0(wgcnaResult, "Beta", Power, "/Net.Rdata"))

  return(net)
}

#' (INTERNAL)
#' @keywords internal

PlotNetwork <- function(EdgeMatrix, NodeIdx, p1, SavePlotFile,
                        EdgeCut = 0.05, VertexLabel = TRUE,
                        AddTrait = TRUE, WeightedEdge = TRUE,
                        ShowNodes = FALSE, PlotTitle = ""){

  # Plot a network of two data types using igraph
  #
  # EdgeMatrix: a matrix indicating the presence and weights (if weighted) of the edges in the network. If a trait is included in the network, it is recorded as the last node for the EdgeMatrix.
  # NodeIdx: a vector indicating node indices.
  # SavePlotFile: the file name for which the plot is saved as.
  # p1: number of features of the first date type.
  # EdgeCut: threshold for edgeweights.
  # VertexLabel: if label the nodes.
  # AddTrait: if include the trait as an additional node.
  # WeightedEdge: whether to use different edge thickness to represent edge weights.
  # ShowNodes: whether to show nodes.
  # PlotTitle: title of the plot.

  nodes <- colnames(EdgeMatrix)
  grp.memb <- nodes[NodeIdx]
  allidx <- matrix(1:nrow(EdgeMatrix), ncol = 1)
  rownames(allidx) <- rownames(EdgeMatrix)
  NodeInfo <- data.frame(id = grp.memb, idx = allidx[grp.memb, ])

  if(AddTrait){
    trait.idx <- nodes[nrow(EdgeMatrix)]
    M.idx <- c(grp.memb, trait.idx)
  }else{
    M.idx <- grp.memb
  }
  M <- as.matrix(EdgeMatrix[M.idx, M.idx])
  M[which(abs(M) < EdgeCut)] <- 0

  net <- graph_from_adjacency_matrix(M, weighted = TRUE,
                                     diag = FALSE, mode = "undirected")
  k <- length(grp.memb)

  vcol <- rep("purple", k); vcol[which(NodeInfo$idx > p1)] <- "orange"
  vshape <- rep("square", k); vshape[which(NodeInfo$idx > p1)] <- "circle"
    if(AddTrait){
      vcol[k+1] <- "SkyBlue"
      vshape[k+1] <- "none"
    }
    lcol <- "blue"
    if(ShowNodes){
      lcol <- "blue"
      if(grp.memb[k] == nodes[nrow(EdgeMatrix)]){
        vcol[k] <- "SkyBlue"
        vshape[k] <- "none"
      }

    }else{
      lcol <- vcol
      vshape <- "none"
    }

    l <- layout_in_circle(net)
    if(!VertexLabel){M.idx <- ""}
    ecol <- rep("gray80", ecount(net))
    if(WeightedEdge){
      ew <- edge.attributes(net)$weight * 5
      ecol[which(ew < 0)] <- "red"
      ew <- abs(ew)
    }else{
      ew <- 5
    }

    pdf(SavePlotFile)
    par(bg = "white")
    plot(net, vertex.color = vcol, vertex.shape = vshape,
         vertex.label.cex = 1.5, layout = l, rescale = FALSE,
         vertex.size = 8, vertex.label = M.idx, edge.width = ew,
         vertex.label.color = lcol, edge.color = ecol, main = PlotTitle)
    dev.off()

}



# Need to first define mRNAmiRNAmodules and p1

#' (INTERNAL)
#' @keywords internal

PlotReducedNetwork <- function(EdgeMatrix, CorrMatrix, mRNAmiRNAmodules,
                               ModuleIdx, AddCorrSign = TRUE, p1 = 3212, SaveFile, CorM = NULL,
                               EdgeCut = 0.05, VertexLabel = TRUE,
                               AddTrait = FALSE, WeightedEdge = FALSE,
                               ShowNodes = TRUE, PlotTitle = ""){

  grp <- mRNAmiRNAmodules[[ModuleIdx]]
  grp.memb <- colnames(EdgeMatrix)[grp]

  p_edge <- nrow(EdgeMatrix)
  p_corr <- nrow(CorrMatrix)

  if(AddTrait){
    if(p_edge != p_corr){
      stop("No edge information for trait. EdgeMatrix and CorrMatrix need to have the same dimension to include a node for trait.")
    }

    trait.id <- colnames(EdgeMatrix)[p_edge]
    M.node <- c(grp.memb, trait.id)

  }else{
    M.node <- grp.memb
  }

  M <- as.matrix(EdgeMatrix[M.node, M.node])
  if(AddCorrSign){M <- M * sign(CorrMatrix[M.node, M.node])}
  M[which(abs(M) < EdgeCut)] <- 0
  newM.node <- M.node[which(apply(abs(M), 1, max) > 0)]

  if(length(newM.node) == 0){
    print("No edge passes threshold.")
  }else{
    M <- M[newM.node, newM.node]
    allidx <- matrix(1:nrow(EdgeMatrix), ncol = 1)
    rownames(allidx) <- rownames(EdgeMatrix)

    NodeInfo <- data.frame(id = newM.node, idx = allidx[newM.node, ])
    net <- graph_from_adjacency_matrix(M, weighted = TRUE,
                                       diag = FALSE, mode = "undirected")

    k <- length(newM.node)
    vcol <- rep("purple", k); vcol[which(NodeInfo$idx > p1)] <- "orange"
    vshape <- rep("square", k); vshape[which(NodeInfo$idx > p1)] <- "circle"
    if(AddTrait){
      vcol[k+1] <- "SkyBlue"
      vshape[k+1] <- "none"
    }
    lcol <- "blue"
    if(ShowNodes){
      lcol <- "blue"
    }else{
      lcol <- vcol
      vshape <- "none"
    }

    if(!VertexLabel){newM.node <- ""}
    ecol <- rep("gray80", ecount(net))
    if(WeightedEdge){
      ew <- abs(edge.attributes(net)$weight) * 5
      ecol[which(edge.attributes(net)$weight < 0)] <- "red"
    }else{
      ew <- 5
    }

    l <- layout_in_circle(net)

    pdf(SaveFile)
    par(bg = "white")
    plot(net, vertex.color = vcol, vertex.shape = vshape, vertex.label.cex = 1.5,
         layout = l, rescale = FALSE, vertex.size = 8, vertex.label = newM.node,
         edge.width = ew, vertex.label.color = lcol, edge.color = ecol)
    title(PlotTitle, cex.main = 1)
    dev.off()
  }

}


############

#' (INTERNAL)
#' @keywords internal

PlotNetwork2 <- function(EdgeMatrix, CorrMatrix, mRNAmiRNAmodules, ModuleIdx,
                         SaveFile, CorM = NULL,
                         EdgeCut = 0.05, VertexLabel = TRUE,
                         AddTrait = FALSE, WeightedEdge = FALSE,
                         ShowNodes = TRUE, PlotTitle = ""){

  # grp <- mRNAmiRNAmodules[ModuleIdx]
  # grp.memb <- colnames(EdgeMatrix)[grps[[mRNAmiRNAmodules[ModuleIdx]]]]
  grp <- mRNAmiRNAmodules[[ModuleIdx]]
  grp.memb <- colnames(EdgeMatrix)[grp]

  allidx <- matrix(1:nrow(EdgeMatrix), ncol = 1)
  rownames(allidx) <- rownames(EdgeMatrix)
  NodeInfo <- data.frame(id = grp.memb, idx = allidx[grp.memb, ])
  if(AddTrait){
    trait.idx <- colnames(EdgeMatrix)[nrow(EdgeMatrix)]
    M.idx <- c(grp.memb, trait.idx)
  }else{
    M.idx <- grp.memb
  }
  M <- as.matrix(EdgeMatrix[M.idx, M.idx])
  M[which(abs(M) < EdgeCut)] <- 0

  M.corr <- CorrMatrix[M.idx, M.idx]
  # M.corr[which(abs(M) < EdgeCut)] <- 0

  net <- graph_from_adjacency_matrix(M, weighted = TRUE,
                                     diag = FALSE, mode = "undirected")
  net.corr <- graph_from_adjacency_matrix(M.corr, weighted = TRUE,
                                          diag = FALSE, mode = "undirected")

  k <- length(grp.memb)

  vcol <- rep("purple", k); vcol[which(NodeInfo$idx > p1)] <- "orange"
  vshape <- rep("square", k); vshape[which(NodeInfo$idx > p1)] <- "circle"
  if(AddTrait){
    vcol[k+1] <- "SkyBlue"
    vshape[k+1] <- "none"
  }
  lcol <- "blue"
  if(ShowNodes){
    lcol <- "blue"
  }else{
    lcol <- vcol
    vshape <- "none"
  }
  # l <- layout.fruchterman.reingold(net)
  # l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
  l <- layout_in_circle(net)
  if(!VertexLabel){M.idx <- ""}
  ecol <- rep("gray80", ecount(net))
  if(WeightedEdge){
    ew <- edge.attributes(net)$weight * 5
    ecol[which(edge.attributes(net.corr)$weight < 0)] <- "red"
  }else{
    ew <- 5
  }

  pdf(SaveFile)
  par(bg = "white")
  plot(net, vertex.color = vcol, vertex.shape = vshape, vertex.label.cex = 1.5,
       layout = l, rescale = FALSE, vertex.size = 8, vertex.label = M.idx,
       edge.width = ew, vertex.label.color = lcol, edge.color = ecol)
  title(PlotTitle, cex.main = 1)
  dev.off()
}


