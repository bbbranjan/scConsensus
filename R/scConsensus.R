#' @title plots the contingency table of two clustering results
#'
#' @author ranjanb
#'
#' @param cluster_labels_1 First vector of cluster labels for each cell
#' @param cluster_labels_2 Second vector of cluster labels for each cell
#' @param automateConsensus Boolean indicating whether automated consensus should be returned
#' @param filename name of contingency table file
#'
#' @return consensus cluster labels vector
plotContingencyTable <- function(cluster_labels_1 = NULL, cluster_labels_2 = NULL, automateConsensus = T, filename = "Contingency_Table.pdf") {

    if(is.null(cluster_labels_1) | is.null(cluster_labels_2)) {
        stop("Incomplete parameters provided.")
    }

    ctg_table <- table(cluster_labels_1, cluster_labels_2)
    ctg_df <- as.data.frame(ctg_table)
    ctg_df <- reshape2::dcast(data = ctg_df, formula = cluster_labels_1~cluster_labels_2, fun.aggregate = sum, value.var = "Freq")
    rownames(ctg_df) <- ctg_df$cluster_labels_1
    ctg_df$cluster_labels_1 <- NULL
    ctg_matrix <- as.matrix(ctg_df)


    pdf(filename, width=15, height=15)

    colorScheme <-
        circlize::colorRamp2(
            seq(-max(abs(
                ctg_matrix
            )), max(abs(
                ctg_matrix
            )), length.out = 5),
            c(
                "cyan",
                "#7FFF7F",
                "yellow",
                "#FF7F00",
                "red"
            )
        )

    ht <- ComplexHeatmap::Heatmap(matrix = ctg_matrix,
                  col = colorScheme,
                  column_title = "Cluster labels 2",
                  row_title = "Cluster Labels 1",
                  column_title_side = "top",
                  column_title_gp = gpar(fontsize = 30, fontface = "bold"),
                  column_names_side = "top",
                  column_names_gp = gpar(fontsize = 20, fontface = "bold"),
                  show_column_dend = FALSE,
                  row_title_side = "left",
                  row_title_gp = gpar(fontsize = 30, fontface = "bold"),
                  row_names_side = "left",
                  row_names_gp = gpar(fontsize = 20, fontface = "bold"),
                  show_row_dend = FALSE,
                  name = "Contingency Table",
                  show_heatmap_legend = FALSE,
                  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                      grid.text(ctg_matrix[i, j], x, y, gp = gpar(fontsize = 20, fontface = "bold", col = "black"))
                  })
    ComplexHeatmap::draw(ht)
    dev.off()

    if(automateConsensus) {
        if(length(cluster_labels_1) > length(cluster_labels_2)){
            consensusClusterLabels = cluster_labels_1
            remainderLabels = cluster_labels_2
        } else if(length(cluster_labels_1) < length(cluster_labels_2)) {
            consensusClusterLabels = cluster_labels_2
            remainderLabels = cluster_labels_1
        } else {
            if(median(table(cluster_labels_1)) > median(table(cluster_labels_2))) {
                consensusClusterLabels = cluster_labels_1
                remainderLabels = cluster_labels_2
            } else {
                consensusClusterLabels = cluster_labels_2
                remainderLabels = cluster_labels_1
            }
        }

        if (ncol(ctg_matrix) > nrow(ctg_matrix)) {
            ctg_matrix <- t(ctg_matrix)
        }

        for(i in 1:nrow(ctg_matrix)) {

            row <- ctg_matrix[i, ]

            row <- 100*(row/sum(row))

            for(j in 1:length(row)) {
                if(overlap >= 10) {
                   consensusClusterLabels[which(remainderLabels) == names(row)[j]] <- paste(rownames(ctg_matrix)[i], names(overlap), sep = "_")
                }
            }

        }
        return(consensusClusterLabels)
    }
}

#' @title Recluster consensus clusters using DE gene analysis
#'
#' @author ranjanb
#'
#' @param dataMatrix the log-transformed and normalized scRNAseq genes x cells matrix
#' @param consensusClusterLabels consensus cell type labels for each cell
#' @param method Method used for DE gene statistical test
#' @param meanScalingFactor scale of the mean gene expression across the gene expression matrix to set a minimum threshold of average cluster expression for a gene to be considered a DE gene.
#' @param qValThrs maximum q-value threshold for statistical test
#' @param fcThrs minimum fold-change threshold for DE gene criterion
#' @param deepSplitValues vector of WGCNA tree cutting deepsplit parameters
#' @param minClusterSize specifies the type of minimum cluster size
#' @param filename name of DE gene object file
#' @param plotName name of DE Heatmap file
#'
#' @return list containing vector of DE genes, clustering tree and dynamic color list.
reclusterDEConsensus <- function(dataMatrix,
                                 consensusClusterLabels,
                                 method = "Wilcoxon",
                                 meanScalingFactor = 5,
                                 qValThrs,
                                 fcThrs,
                                 deepSplitValues = 1:4,
                                 minClusterSize = 10,
                                 filename = "de_gene_object.rds",
                                 plotName = "DE_Heatmap") {

    ### Initialise data input
    dataIn <- as.matrix(dataMatrix)


    ### Initialize mean expression threshold
    meanExprsThrs = meanScalingFactor * mean(expm1(dataIn))

    ### Use only clusters with number of cells > minimum cluster size for DE gene calling
    which(table(consensusClusterLabels) > minClusterSize)
    colorCounts <-
        table(consensusClusterLabels)

    ### Extract unique cluster labels
    uniqueClusters <-
        names(colorCounts[colorCounts > minClusterSize])

    # ignore the grey cluster because it represents unclustered cells
    uniqueClusters <-
        uniqueClusters[!grepl("grey", uniqueClusters)]

    ### If the cluster label vector is unnamed, name it in the order of the data matrix columns
    if (is.null(names(consensusClusterLabels))) {
        names(consensusClusterLabels) <- colnames(dataIn)
    }

    ### Initialize nested lists to store q values, log-normalized fold change values and de genes
    qValueList <-
        rep(list(list()), length(uniqueClusters))
    log2FCList <-
        rep(list(list()), length(uniqueClusters))
    deGeneList <- rep(list(list()), length(uniqueClusters))

    ### Initialise number of comparisons as n(n-1)/2
    numComparisons <-
        (length(uniqueClusters) * (length(uniqueClusters) - 1)) / 2

    ### Conduct pairwise cluster comparison to obtain q-values, log-normalized fold changes and
    ### DE genes for each comparison
    for (i in 1:(length(uniqueClusters) - 1)) {
        for (j in (i + 1):length(uniqueClusters)) {
            ### Get the cell names cell data for cluster i
            cellNamesi <-
                names(consensusClusterLabels)[which(consensusClusterLabels == uniqueClusters[i])]
            cellDatai <- dataIn[, cellNamesi]

            ### Get the cell names cell data for cluster j
            cellNamesj <-
                names(consensusClusterLabels)[which(consensusClusterLabels == uniqueClusters[j])]
            cellDataj <- dataIn[, cellNamesj]

            deCellData <- cbind(cellDatai, cellDataj)

            ### Declare variables to store p-values, q-values and log2 fold change values for each pairwise gene comparison
            pval <- NA
            qval <- NA
            log2fc <- NA

            if (method == "Wilcoxon") {
                ### For each gene, conduct Wilcoxon Rank Sum test to obtain q-value and log2-normalized fold change
                for (k in 1:nrow(cellDatai)) {

                    # extract data of gene k for cluster i
                    geneDatai <- cellDatai[k, ]

                    # extract data of gene k for cluster j
                    geneDataj <- cellDataj[k, ]

                    # run wilcoxon test
                    wcTestOut <-
                        wilcox.test(geneDatai, geneDataj)

                    # initialize q value and fold change for each gene
                    pval[k] <- wcTestOut$p.value
                    log2fc[k] <-
                        mean(geneDatai) - mean(geneDataj)
                }

                # check if gene satisfies mean expression threshold
                meanExprsLogicalVector <-
                    apply(deCellData, 1, function(row) {
                        mean(row[cellNamesi]) > log2(meanExprsThrs) |
                            mean(row[cellNamesj]) > log2(meanExprsThrs)
                    })

                # adjust p-value using Benjamini-Hochberg adjustment
                qval <-
                    p.adjust(
                        p = pval,
                        method = "BH",
                        n = nrow(cellDatai)
                    )

            } else if (method == "edgeR") {

                # check if gene satisfies mean expression threshold
                meanExprsLogicalVector <-
                    apply(deCellData, 1, function(row) {
                        mean(row[cellNamesi]) > log2(meanExprsThrs) |
                            mean(row[cellNamesj]) > log2(meanExprsThrs)
                    })

                # Create DGE object
                dgeObj <- edgeR::DGEList(counts=deCellData, group = c(rep(1, ncol(cellDatai)), rep(-1, ncol(cellDataj))))

                # Estimate common and tag-wise dispersion for the pair of clusters
                dgeObj <- edgeR::estimateCommonDisp(dgeObj)
                dgeObj <- edgeR::estimateTagwiseDisp(dgeObj)

                # Calculate norm factors for data
                dgeObj <- edgeR::calcNormFactors(dgeObj, method = "none")

                # Perform exact test
                etObj <- edgeR::exactTest(object = dgeObj)

                # Save p-value results
                pval <- etObj$table$PValue

                # adjust p-value using Benjamini-Hochberg adjustment
                qval <-
                    p.adjust(
                        p = pval,
                        method = "BH",
                        n = nrow(cellDatai)
                    )

                log2fc <- etObj$table$logFC

            } else {
                ### exit from function if no method is chosen
                print("Incorrect method chosen.")
                return(NULL)
            }

            ### determine if a gene is a DE gene based on thresholds
            deGeneLogicalVector <-
                qval < qValThrs &
                abs(log2fc) > log2(fcThrs)

            # filter de genes by mean expression
            deGeneLogicalVector <- deGeneLogicalVector & meanExprsLogicalVector

            print(paste0(
                uniqueClusters[i],
                ", ",
                uniqueClusters[j],
                " DE genes: ",
                sum(deGeneLogicalVector)
            ))

            ### store q-values, log-normalized fold change and DE genes
            qValueList[[i]][[j]] <- qval
            log2FCList[[i]][[j]] <- log2fc
            deGeneList[[i]][[j]] <-
                rownames(cellDatai)[deGeneLogicalVector]

            if (sum(deGeneLogicalVector) <= 1)
                next

            de_genes <- deGeneList[[i]][[j]]


        }
    }

    ### Name the q-value, log2 fold change and DE gene lists by cluster names
    names(qValueList) <- uniqueClusters
    names(log2FCList) <- uniqueClusters
    names(deGeneList) <- uniqueClusters

    saveRDS(object = qValueList, file = "qValueList.rds")
    saveRDS(object = log2FCList, file = "log2FCList.rds")
    saveRDS(object = deGeneList, file = "deGeneList.rds")

    ### Obtain union of all de genes
    # initialize empty DE gene union vector
    deGeneUnion <- c()

    # For each pair of clusters
    for (i in 1:(length(uniqueClusters) - 1)) {
        for (j in (i + 1):length(uniqueClusters)) {

            # Rank the DE genes for each comparison by fold change and take the top 30

            de_log2fc <- log2FCList[[i]][[j]][which(rownames(dataIn) %in% deGeneList[[i]][[j]])]
            names(de_log2fc) <- rownames(dataIn)[which(rownames(dataIn) %in% deGeneList[[i]][[j]])]
            sorted_de_log2fc <- sort(x = abs(de_log2fc), decreasing = T)

            if(length(sorted_de_log2fc) > 30) {
                de_genes <- names(sorted_de_log2fc)[1:30]
        } else {
                de_genes <- names(sorted_de_log2fc)
            }

            deGeneUnion <-
                union(deGeneUnion, de_genes)
        }
    }

    print(str(deGeneUnion))

    saveRDS(deGeneUnion, file = "deGeneUnion.rds")

    ### Compute PCA + euclidean distance of cells based on the union of DE genes
    pca.data <- irlba::prcomp_irlba(x = t(dataIn[deGeneUnion, ]), n = min(length(deGeneUnion), 15), center = TRUE, scale. = FALSE)$x

    d <- dist(pca.data, method = "euclidean")

    ### Compute 2-level correlation distance of cells based on the union of DE genes
    # d = as.dist(1 - cor(dataIn[deGeneUnion, ], method = "pearson"))

    ### Build dendrogram using distance matrix d
    cellTree = stats::hclust(d, method = "ward.D2")

    # Initialize empty dynamic colors list
    dynamicColorsList <- list()

    ### For all values of deepsplit specified
    for (dsv in deepSplitValues) {
        ### Use WGCNA dynamic tree cutting to obtain clusters
        dynamicGroups = dynamicTreeCut::cutreeDynamic(
            dendro = cellTree,
            distM = as.matrix(d),
            deepSplit = dsv,
            pamStage = FALSE,
            minClusterSize = minClusterSize
        )
        dynamicColors = WGCNA::labels2colors(dynamicGroups)

        dynamicColorsList[[paste("deepsplit:", dsv)]] <-
            dynamicColors

    }

    ### Name dynamic colors list
    names(dynamicColorsList) <- paste("deepsplit:", deepSplitValues)

    ### Compute number of detected genes for each cell in the input data
    nodg <- c()
    for (i in 1:length(colnames(dataIn))) {
        nodg[i] <- length(which(dataIn[, i] > 0))
    }

    ### Create and initialise object to return to function caller
    returnObj = list(
        "deGeneUnion" = deGeneUnion,
        "cellTree" = cellTree,
        "dynamicColors" = dynamicColorsList
    )

    ### Save DE object
    saveRDS(object = returnObj, file = filename)

    ### Plot DE Gene Plot
    cellTypeDEPlot(
        dataMatrix = dataIn[deGeneUnion,],
        nodg = nodg,
        cellTree = cellTree,
        clusterLabels = consensusClusterLabels,
        dynamicColorsList = dynamicColorsList,
        colScheme = "violet",
        filename = plotName
    )

    return(returnObj)
}



