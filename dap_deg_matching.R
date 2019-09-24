



run_sigma <- function(rna_sig, atac_sig, range_bp = c(1e4, 2e4, 1e5, 2e5, 5e5), gene_id = "name", peak_id = "orig_name") {
    peak_pool_up <- sapply(atac_sig, function(x) {
        length(x)
    })
    
    peak_all  <- unique(unlist(as(atac_sig, "GRangesList")))
    
    gene_pool_up <- sapply(rna_sig, function(x) {
        #length(x)
        length(unique(mcols(x)[[gene_id]]))
    })
    gene_all <- unique(unlist(as(rna_sig, "GRangesList")))
    
    pg_list <- list()
    hyp_deg_list <- list()
    hyp_dap_list <- list()
    for(range in range_bp) {
        print(range)
        
        # DAP enrichment near DEG
        pg_overlap_list<- lapply(atac_sig, function(p) {
            lapply(rna_sig, function(t){
                findOverlaps(p, t, maxgap = range)
                pt_hit <- findOverlaps(p, t, maxgap = range)
                hit_tbl <- as.data.frame(pt_hit)
                hit_tbl$DAP <- mcols(p[hit_tbl$queryHits])[[peak_id]]
                hit_tbl$DEG <- mcols(t[hit_tbl$subjectHits])[[gene_id]]
                return(hit_tbl)
            })
        })
        
        overlap_length_peak<-as.data.frame(lapply(pg_overlap_list, function(x) {
            sapply(x, function(y){
                length(unique(y$DAP))
            })
        }))
        
        peak_all_og <- sapply(rna_sig, function(t){
            peak_og<-findOverlaps(peak_all, t, maxgap = range)
            hit_tbl <- as.data.frame(peak_og)
            hit_tbl$DAP <- mcols(peak_all[hit_tbl$queryHits])[[peak_id]]
            hit_tbl$DEG <- mcols(t[hit_tbl$subjectHits])[[gene_id]]
            length(unique(hit_tbl$DAP))
        })
        
        peak_dif_og <- sapply(rna_sig, function(t){
            peak_og<-findOverlaps(peak_all, t, maxgap = range)
            hit_tbl <- as.data.frame(peak_og)
            hit_tbl$DAP <- mcols(peak_all[hit_tbl$queryHits])[[peak_id]]
            hit_tbl$DEG <- mcols(t[hit_tbl$subjectHits])[[gene_id]]
            length(unique(mcols(peak_all)[[peak_id]])) - length(unique(hit_tbl$DAP))
        })
        
        hyp_dap<-sapply(1:ncol(overlap_length_peak), function(i) {
            x <- overlap_length_peak[,i]
            sapply(1:length(x), function(j){
                phyper(x[j], peak_all_og[j], peak_dif_og[j], peak_pool_up[i], lower.tail = F)
            })
        })
        colnames(hyp_dap) <- colnames(overlap_length_peak)
        rownames(hyp_dap) <- rownames(overlap_length_peak)
        hyp_dap_list[[as.character(range)]] <- hyp_dap
        pg_list[[as.character(range)]] <- pg_overlap_list
        
        
        # DEG enrichment near DAP
        overlap_length_gene<-as.data.frame(lapply(pg_overlap_list, function(x) {
            sapply(x, function(y){
                length(unique(y$DEG))
            })
        }))
        gene_all_op <- sapply(atac_sig, function(p){
            gene_og<-findOverlaps(gene_all, p, maxgap = range)
            hit_tbl <- as.data.frame(gene_og)
            hit_tbl$DEG <- mcols(gene_all[hit_tbl$queryHits])[[gene_id]]
            hit_tbl$DAP <- mcols(p[hit_tbl$subjectHits])[[peak_id]]
            length(unique(hit_tbl$DEG))
        })
        gene_dif_op <- sapply(atac_sig, function(p){
            gene_og<-findOverlaps(gene_all, p, maxgap = range)
            hit_tbl <- as.data.frame(gene_og)
            hit_tbl$DEG <- mcols(gene_all[hit_tbl$queryHits])[[gene_id]]
            hit_tbl$DAP <- mcols(p[hit_tbl$subjectHits])[[peak_id]]
            length(unique(mcols(gene_all)[[gene_id]])) - length(unique(hit_tbl$DEG))
        })
        hyp_deg<-sapply(1:nrow(overlap_length_gene), function(i) {
            x <- overlap_length_gene[i,]
            sapply(1:length(x), function(j){
                phyper(unlist(x[j]), gene_all_op[j], gene_dif_op[j], gene_pool_up[i], lower.tail = F)
            })
        })
        hyp_deg <- t(hyp_deg)
        colnames(hyp_deg) <- colnames(overlap_length_gene)
        rownames(hyp_deg) <- rownames(overlap_length_gene)
        hyp_deg_list[[as.character(range)]] <- hyp_deg
    }
    return(list(dap = hyp_dap_list, deg = hyp_deg_list, tbls = pg_list))
}



plot_matching <- function(matching_res, pval_min = 1e-50, scale = "none", plotFile = NULL, width, height) {
    if(scale == "x") {
        scale_func <- function(x) t(scale(t(x)))
    } else if(scale == "y") {
        scale_func <- function(x) scale(x)
    } else {
        scale_func <- function(x) x
    }
    dap_hmaps<-lapply(1:length(matching_res[["dap"]]), function(i) {
        x <- matching_res[["dap"]][[i]]
        range <- names(matching_res[["dap"]])[i]
        pres<-pheatmap(scale_func(-log10(x+pval_min)), cluster_rows = F, cluster_cols = F, main = paste0("DAP_wise_distance_range_", range))
        pres[[4]]
    })
    deg_hmaps<-lapply(1:length(matching_res[["deg"]]), function(i) {
        x <- matching_res[["deg"]][[i]]
        range <- names(matching_res[["deg"]])[i]
        pres<-pheatmap(scale_func(-log10(x+pval_min)), cluster_rows = F, cluster_cols = F, main = paste0("DEG_wise_distance_range_", range))
        pres[[4]]
    })
    pdf(plotFile, width = width, height = height)
    do.call(grid.arrange, c(dap_hmaps, ncol=length(dap_hmaps)))
    do.call(grid.arrange, c(deg_hmaps, ncol=length(deg_hmaps)))
    dev.off()
    dev.off()
}



















