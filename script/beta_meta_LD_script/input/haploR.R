
if ( require("haploR") ) {
    print("TRUE: The haploR library exists")
} else {
    print("FALSE: The haploR library doesn't exist")
    install.packages("haploR", repos = "http://cran.us.r-project.org")
}


if ( require("tidyverse") ) {
    print("TRUE: The tidyverse library exists")
} else {
    print("FALSE: The tidyverse library doesn't exist")
    install.packages("tidyverse")
}


library(tidyverse)
getCurrentFileLocation <-  function()
{
    this_file <- commandArgs() %>% 
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
    if (length(this_file)==0)
    {
      this_file <- rstudioapi::getSourceEditorContext()$path
    }
    return(dirname(this_file))
}

input_txt <- paste0(getCurrentFileLocation() ,"/input_SNPs.txt")
rsIDs <- readLines(input_txt)

for (x in rsIDs) {
  xx <- queryHaploreg(query=x, ldThresh=1, ldPop="ASN")
  # print(xx)
  subset.LD1 <- xx[c("rsID")]
  # print(subset.LD1)
  txtfile <- paste0(str_sub(getCurrentFileLocation(), start = 1, end = -7), "/output/correlated_with_",x,".txt")
  write.table(subset.LD1, file=txtfile, quote=FALSE, sep="\t")
}
