# CONSTANTS

# data files directory 
# make sure these directories exist
FDIR = "../data_files"
OUTDIR = "../out_files"

# separation for file names
FSEP = "/"

# Wald test significance value
ALPHA = 0.05

compare_dats <- function(dat1, dat2){

    #' Calculate the difference matrix

    # calculates difference matrices
    diff = dat1 - dat2
    # diff_err = sqrt(err1^2 + err2^2)

    return (diff)
    # return diff_err
}

calc_wald <- function(dmat, err, alpha){

    #'Calculate matrix of Wald values given values and errors

    # null hypothesis
    p = 0

    # if no error matrix given
    # W <- (dmat-p)^2 / (dmat*(1-dmat)/n) 

    # calculate score
    W <- (dmat-p)^2 / err
    pval = 1-mapply(function (x) pchisq(x, df=1), W)    # p-value 

    return (W)

}

compress_to_pos <- function(mat){

    #' Compress to array at position

    return (colMeans(mat, na.rm=TRUE))
}

to_csv <- function(fname){

    #' Turns filename string into filename.csv

    return (paste(fname, 'csv', sep='.'))
}

run <- function(fname1, fname2, err){

    #' Main run 
    fname1 = "test1"
    fname2 = "test2"
    err="testerr"

    # load data (no row or column names)
    dat1 = data.matrix(read.csv(paste(FDIR, to_csv(fname1), sep=FSEP), header=TRUE, check.names=FALSE), rownames.force = NA) # row.names=1
    dat2 = data.matrix(read.csv(paste(FDIR, to_csv(fname2), sep=FSEP), header=TRUE, check.names=TRUE), rownames.force = NA)

    err_mat = data.matrix(read.csv(paste(FDIR, to_csv(err), sep=FSEP), header=TRUE, check.names=FALSE), rownames.force = NA)

    # get difference matrix
    diff_mat = compare_dats(dat1, dat2)

    # get diff by position
    diff_lst = compress_to_pos(diff_mat)

    # do stat tests
    wald_mat = calc_wald(diff_mat, err_mat, ALPHA)

    # get wald scores by pos
    wald_lst = compress_to_pos(wald_mat)

    # write results to matrix
    write.table(diff_mat, file = paste(OUTDIR, to_csv(paste(fname1, fname2, 'diff_mat', sep='_')), sep="/"), sep=",")
    write.table(diff_lst, file = paste(OUTDIR, to_csv(paste(fname1, fname2, 'diff_lst', sep='_')), sep="/"), sep=",")
    write.table(wald_mat, file = paste(OUTDIR, to_csv(paste(fname1, fname2, 'wald_mat', sep='_')), sep="/"), sep=",")
    write.table(wald_lst, file = paste(OUTDIR, to_csv(paste(fname1, fname2, 'wald_lst', sep='_')), sep="/"), sep=",")


    # ========================
    # IF PLOTTING, UNCOMMENT
    # ========================
    # heatmaps of differences
    # change cluster_cols = TRUE to cluster by position

    # library("pheatmap")
    # pheatmap(diff_mat, cluster_rows=FALSE, cluster_cols=FALSE, main="Difference Heatmap", show_colnames=TRUE )

    # histogram of difference by position

    # hist(diff_lst)

    # heatmaps of alpha values

    # library("pheatmap")
    # pheatmap(wald_mat, cluster_rows=FALSE, cluster_cols=FALSE, main="Difference Heatmap", show_colnames=TRUE )

    # histogram of wald test score by position

    # hist(wald_lst)


}

# ========================
# test run
# ========================

run("test1", "test2", "testerr")

# ========================
# actual run
# ========================

# run("fitness1", "fitness2", "err")



