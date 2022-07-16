#' Read a vcf file
#'
#' @param vcfFile Path to vcf file
#'
#' @return A list with elements header, vm, vi
#' @export
#'
#' @examples
#' \dontrun{
#' vcf<-read_vcf("inst/extdata/trio.2010_06.ychr.sites.vcf.gz")
#' }

read_vcf<-function(vcfFile) {

    vv=readr::read_tsv(vcfFile,comment="##",show_col_types = FALSE) %>%
        rename(CHROM=`#CHROM`) %>%
        mutate(VID=row_number()) %>%
        select(VID,everything())

    vs=select(vv,VID,10:ncol(vv)) %>%
        group_split(FORMAT) %>%
        purrr::map(parse_vcf_samples) %>%
        bind_rows %>%
        arrange(VID,SAMPLE)

    vm=select(vv,1:9)

    v0=vm %>%
        select(VID,INFO) %>%
        separate_rows(INFO,sep=";") %>%
        filter(!grepl("=",INFO)) %>% mutate(VALUE="T")

    v0flds=distinct(v0,INFO) %>% pull

    v1=vm %>%
        select(VID,INFO) %>%
        separate_rows(INFO,sep=";") %>%
        filter(grepl("=",INFO)) %>%
        separate(INFO,c("INFO","VALUE"),sep="=")

    vi=bind_rows(v0,v1) %>% spread(INFO,VALUE) %>% mutate_at(v0flds,~ifelse(is.na(.),FALSE,TRUE))
    vm=vm %>% select(-INFO) %>% full_join(vi,by="VID")

    ni=2
    nj=3

    while(!any(grepl("^#CHROM",readLines(vcfFile,n=nj)))) {
        ii=nj
        nj=ni+nj
        ni=ii
    }

    header=readLines(vcfFile,n=nj)
    hi=grep("^#CHROM",header)
    header=header[1:(hi-1)]

    list(header=header,vm=vm,vs=vs)

}

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
})

parse_vcf_samples<-function(vi) {

    ff=strsplit(vi$FORMAT[1],":")[[1]]
    vi %>%
        select(-FORMAT) %>%
        gather(SAMPLE,GDATA,-VID) %>%
        separate(GDATA,ff,sep=":",fill="right") %>%
        arrange(VID)

}
