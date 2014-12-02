mat <- matrix(rexp(200, rate=.1), ncol=20)

write.table(mat, file="testerr.csv", sep=",")