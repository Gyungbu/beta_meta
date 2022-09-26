if ( require("haploR") ) {
    print("TRUE: The library exists")
} else {
    print("FALSE: The library doesn't exist")
    install.packages("haploR", repos = "http://cran.us.r-project.org")
}

input_txt <- paste0(getwd(),"/input/input_SNPs.txt")
rsIDs <- readLines(input_txt)

for (x in rsIDs) {
  xx <- queryHaploreg(query=x, ldThresh=1, ldPop="ASN")
  # print(xx)
  subset.LD1 <- xx[c("rsID")]
  # print(subset.LD1)
  txtfile <- paste0(getwd(),"/output/correlated_with_",x,".txt")
  write.table(subset.LD1, file=txtfile, quote=FALSE, sep="\t")
}