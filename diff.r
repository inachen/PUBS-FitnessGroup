compare_dats <- function(dat1, dat2, err1, err2){
    # calculates difference and graphs heatmap

    diff = dat1 - dat2
    diff_err = sqrt(err1^2 + err2^2)

    library("pheatmap")

    # heatmaps of differences
    # change cluster_cols = TRUE to cluster by position
    pheatmap(diff, cluster_rows=FALSE, cluster_cols=FALSE, main="Difference Heatmap", show_colnames=TRUE )
    pheatmap(diff_err, cluster_rows=FALSE, cluster_cols=FALSE, main="Error Heatmap", show_colnames=TRUE )
}

calc_scores <- function(dat1, dat2, err1, err2){

    # TO-DO!
    # check out t-tests and what not that r have functions for

}

run <- function(){
    # main run function

    # load data
    # change path and parameters as needed
    dat1 = read.csv("data_files/fitness1.csv", row.names=1, header=TRUE, check.names=FALSE)
    dat2 = read.csv("data_files/fitness2.csv", row.names=1, header=TRUE, check.names=FALSE)

    # if err values given in different files
    # may just need to extract err values from the data csv's
    err1 = read.csv("data_files/err1.csv", row.names=1, header=TRUE, check.names=FALSE)
    err2 = read.csv("data_files/err2.csv", row.names=1, header=TRUE, check.names=FALSE)

    # draw heatmaps
    compare_dats(dat1, dat2, err1, err2)

    # do stat tests
    calc_scores(dat1, dat2, err1, err2)
}