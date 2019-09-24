

library(clusterProfiler)
library(reldist)

ifg_select <- function(data, cluster, cluster_min_cell_num = 100, min_cluster_expr_fraction = .1, gini_cut = NULL, gini_cut_qt = .75, compute_go = F, filePath = NULL, fileName = "ifg_select_", orgdb = "org.Mm.eg.db", gene_id_type = "ENSEMBL", gene_background = NULL, return_all = F) {
    if(ncol(data) != length(cluster)) stop("Cell number do not match cluster length.")
    use_clus <- names(which(table(cluster) >= cluster_min_cell_num))
    # For each cluster, compute gene expressed fraction
    expr_clus_frac<-sapply(use_clus, function(x) {
        cur_data <- data[,cluster == x]
        Matrix::rowMeans(cur_data > 0)
    })
    # Compute gini coefficient 
    # Require a gene to be expressed in at least one cluster with at least .1 expressed fraction to be considered for downstream uniform gene selection
    use_g <- rownames(expr_clus_frac)[rowSums(expr_clus_frac >= min_cluster_expr_fraction) > 0] # 11813
    message(paste0("Selecting informative features from ", length(use_g), " robustly detected features."))
    expr_clus_frac <- expr_clus_frac[rownames(expr_clus_frac) %in% use_g,]

    gene_clus_gini <- apply(expr_clus_frac, 1, gini)

    if(is.null(gini_cut)) {
        gini_cut <- quantile(gene_clus_gini, gini_cut_qt) # 0.394804 
        message(paste0("Cut at gini quantile ", gini_cut_qt, " with value ", gini_cut))
    } else {
        message(paste0("Cut at gini value ", gini_cut))
    }
    
    pdf(paste0(filePath, fileName, "gene_clus_gini_hist.pdf"))
    hist(gene_clus_gini, breaks = 100)
    abline(v = gini_cut, col=c("red"), lty=c(2), lwd=c(3))
    dev.off()
    exclude_g <- names(gene_clus_gini)[gene_clus_gini < gini_cut]
    include_g <- names(gene_clus_gini)[gene_clus_gini >= gini_cut]
    write.csv(data.frame(less_specific_feature = exclude_g), paste0(filePath, fileName, "less_specific_feature_list.csv"))
    write.csv(data.frame(specific_feature = include_g), paste0(filePath, fileName, "specific_feature_list.csv"))
    
    message(paste0("Found ", length(exclude_g), " less specific features."))
    message(paste0("Returning ", length(include_g), " specific features."))
    
    if(compute_go) {
        if(!length(gene_background)) {
            gene_background <- use_g
        }
        exclude_g.df <- bitr(exclude_g, fromType = gene_id_type,
                             toType = c("SYMBOL", "ENTREZID"),
                             OrgDb = orgdb)
        bg.df <- bitr(gene_background, fromType = gene_id_type,
                      toType = c("SYMBOL", "ENTREZID"),
                      OrgDb = orgdb)
        exclude_g_go<- enrichGO(gene        = exclude_g.df$ENTREZID,
                                universe      = bg.df$ENTREZID,
                                OrgDb         = orgdb,
                                ont           = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.05,
                                readable      = TRUE)
        write.csv(exclude_g_go, paste0(filePath, fileName, "less_specific_feature_go.csv"))
        
        
        include_g.df <- bitr(include_g, fromType = gene_id_type,
                             toType = c("SYMBOL", "ENTREZID"),
                             OrgDb = orgdb)
        bg.df <- bitr(gene_background, fromType = gene_id_type,
                      toType = c("SYMBOL", "ENTREZID"),
                      OrgDb = orgdb)
        include_g_go<- enrichGO(gene        = include_g.df$ENTREZID,
                                universe      = bg.df$ENTREZID,
                                OrgDb         = orgdb,
                                ont           = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.05,
                                readable      = TRUE)
        write.csv(include_g_go, paste0(filePath, fileName, "specific_feature_go.csv"))
    }
    if(return_all) {
        return(list(include_g = include_g, exclude_g = exclude_g))
    } else {
        return(include_g)   
    }
}






