# CONSTANTS

# data files directory 
# make sure these directories exist
FDIR = "data_files"
OUTDIR = "out_files"

# separation for file names
FSEP = "/"

# Wald test significance value
ALPHA = 0.05

compare_dats <- function(dat1, dat2){

    '''Calculate the difference matrix'''

    # calculates difference matrices
    diff = dat1 - dat2
    # diff_err = sqrt(err1^2 + err2^2)

    return diff
    # return diff_err
}

calc_wald <- function(dmat, err, alpha){

    '''Calculate matrix of Wald values given values and errors'''

    # null hypothesis
    p = 0

    # if no error matrix given
    # W <- (dmat-p)^2 / (dmat*(1-dmat)/n) 

    # calculate score
    W <- (dmat-p)^2 / err
    pval = 1-pchisq(W, df=1)    # p-value 

}

compress_to_pos <- function(mat){

    '''compress to array at position'''

    return colMeans(mat, na.rm=TRUE)
}

to_csv <- function(fname){

    return paste(fname, 'csv', sep='.')
}

run <- function(fname1, fname2, err){

    # load data (no row or column names)
    dat1 = read.csv(paste(FDIR, to_csv(fname1), sep=FSEP), header=FALSE, check.names=FALSE) # row.names=1
    dat2 = read.csv(paste(FDIR, to_csv(fname2), sep=FSEP), header=FALSE, check.names=FALSE)

    err_mat = read.csv(paste(FDIR, err, sep=FSEP), header=FALSE, check.names=FALSE)

    # get difference matrix
    diff_mat = compare_dats(dat1, dat2)

    # get diff by position
    diff_lst = compress_to_pos(diff_mat)

    # do stat tests
    wald_mat = calc_wald(diff_mat, err_mat, ALPHA)

    # get wald scores by pos
    wald_lst = compress_to_pos(wald_mat)

    # write results to matrix
    write.matrix(diff_mat, file = paste(OUTDIR, to_csv(paste(fname1, fname2, 'diff_mat', '_')), sep="/"), sep=",")
    write.matrix(diff_lst, file = paste(OUTDIR, to_csv(paste(fname1, fname2, 'diff_lst', '_')), sep="/"), sep=",")
    write.matrix(wald_mat, file = paste(OUTDIR, to_csv(paste(fname1, fname2, 'wald_mat', '_')), sep="/"), sep=",")
    write.matrix(wald_lst, file = paste(OUTDIR, to_csv(paste(fname1, fname2, 'wald_lst', '_')), sep="/"), sep=",")


    #########################
    # IF PLOTTING, UNCOMMENT
    #########################
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

run("fitness1", "fitness2", "err")



