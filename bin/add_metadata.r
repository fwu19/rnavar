#!/usr/bin/env Rscript

# Author: @fwu19

options(stringsAsFactors = F)
options(scipen = 99)
library(dplyr)

## functions ####
## add default value if a column is missing
fill_column <- function(df, colv, default.value, na.value = NULL, missing.value = NULL){
    if (!colv %in% colnames(df)){
        df[,colv] <- default.value
    }
    
    if (!is.null(na.value)){
        df[,colv] <- ifelse(is.na(df[,colv]), na.value, df[,colv])
    }
    
    if (!is.null(missing.value)){
        df[,colv] <- ifelse(df[,colv] == "", missing.value, df[,colv])
    }
    
    return(df)
}

add_metadata <- function(ss, meta_csv){
    ## update with metadata if provided

    if (file_test('-f', meta_csv) & grepl('.csv$', meta_csv) & file.size(meta_csv) > 0){
        meta <- read.csv(meta_csv)
        if (!'sample' %in% colnames(meta)){
            stop ( 'Missing column sample!' )
        }
        
        ss <- ss %>% 
            inner_join(
                meta, by = 'sample', suffix = c(".x", "")
            ) %>% 
            dplyr::select(!ends_with(".x")) 
    }
    
    return(ss)
}

## read arguments ####
args <- as.vector(commandArgs(T))
in_csv <- args[1]
meta_csv <- ifelse(length(args) > 1, args[2], '')
out_csv <- ifelse(length(args) > 2, args[3], in_csv)

## generate sample sheet ####
ss <- read.csv(in_csv)

## check single end fastq
if (!'sample' %in% colnames(ss)){
    stop ( 'Missing column sample!' )
}

## add metadata ####
ss <- add_metadata(ss, meta_csv)
ss %>% 
    write.table(out_csv, sep = ',', quote = F, row.names = F)

