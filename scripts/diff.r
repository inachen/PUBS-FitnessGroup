compare_dats <- function(dat1, dat2, err1, err2){
    # calculates difference and graphs heatmap

    diff = dat1 - dat2
    diff_err = sqrt(err1^2 + err2^2)

    # library("pheatmap")

    # heatmaps of differences
    # change cluster_cols = TRUE to cluster by position
    # pheatmap(diff, cluster_rows=FALSE, cluster_cols=FALSE, main="Difference Heatmap", show_colnames=TRUE )
    # pheatmap(diff_err, cluster_rows=FALSE, cluster_cols=FALSE, main="Error Heatmap", show_colnames=TRUE )


    return diff
}

calc_wald <- function(dmat, err, alpha){

    # calculate wald scores for two matrices 
    p = 0

    # W <- (dmat-p)^2 / (dmat*(1-dmat)/n)
    W <- (dmat-p)^2 / err
    pval = 1-pchisq(W, df=1)    # p-value 

}

compress_to_pos <- function(mat){

    return colMeans(mat, na.rm=TRUE)
}

run <- function(mat){
    # main run function

    # load data
    # change path and parameters as needed
    dat1 = read.csv("data_files/fitness1.csv", row.names=1, header=TRUE, check.names=FALSE)
    dat2 = read.csv("data_files/fitness2.csv", row.names=1, header=TRUE, check.names=FALSE)

    # if err values given in different files
    # may just need to extract err values from the data csv's
    err_mat = read.csv("data_files/err.csv", row.names=1, header=TRUE, check.names=FALSE)

    # draw heatmaps
    diff_mat = compare_dats(dat1, dat2, err1, err2)

    # get diff by position
    diff_lst = compress_to_pos(diff_mat)

    # do stat tests
    alpha = 0.05
    wald_mat = calc_wald(diff_mat, err_mat, alpha)

    # get wald scores by pos
    wald_lst = compress_to_pos(wald_mat)
}

run()



