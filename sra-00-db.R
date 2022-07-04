rm(list = ls())
#' https://www.ncbi.nlm.nih.gov/Traces/wgs/?page=1&view=all
#' 02 July 2022
sra_data <- readr::read_csv(
  '../../data/sra/wgs_selector.csv'
)
sra_data <- sra_data[!is.na(sra_data$bioproject_s),]
sra_data[sra_data$bioproject_s == 'LSAR00000000',]

manual_sra <- readr::read_csv('../../data/manual-isolates/batch-01-antibiotics.csv')

sra_projects_manual <- gsub(
  ' ',
  '',unique(gsub('\\.','',na.omit(manual_sra$bioproject)))
  )
sra_proj_manual <- as.character(unlist(lapply(sra_projects_manual,function(x){
  strsplit(x,split = ',')[[1]]
})))

sra_proj_manual <- as.character(sapply(sra_proj_manual,function(x){
  xx <- strsplit(x,split = '')[[1]]
  paste0(xx[grepl(pattern = "[[:alnum:]]",xx)],collapse = '')
}))

idx_manual <- sra_data$bioproject_s %in% sra_proj_manual

sra_manual <- sra_data[idx_manual,]
setdiff(sra_proj_manual,unique(sra_manual$bioproject_s))
rentrez::entrez_summary(db="bioproject",id=306708)

#---- GET EntezID
eSearch <- list()
for(i in 1:length(sra_proj_manual)){
  eSearch[[i]] <- rentrez::entrez_search(
    db="bioproject",term=sra_proj_manual[i]
  )  
}
eSummary <- list()
for(i in 1:length(eSearch)){
  eSummary[[i]] <- rentrez::entrez_summary(db="bioproject",id=eSearch[[i]]$ids)
}

eProject <- list()
for(i in 1:length(eSearch)){
  print(i)
  xx <- XML::xmlToList(rentrez::entrez_fetch(
    db = "bioproject",id = eSearch[[i]]$ids,rettype = 'xml',
    parsed = TRUE
  ))
  eProject[[i]] <- xx
}
extract_details <- function(x){
  xx_items <- x[[1]]$Project$ProjectDescr
  xx_attr <- names(xx_items)
  locus_tags <- which(xx_attr == "LocusTagPrefix")
  biosamps <- c()
  j <- 1
  for(i in locus_tags){
    biosamps[j] <- xx_items[[i]]$.attrs
    j <- j+1
  }
  xx_df <- data.frame(
    bioproject = x[[1]]$Project$ProjectID$ArchiveID['accession'],
    id = x[[1]]$Project$ProjectID$ArchiveID['id'],
    biosamples = biosamps,
    db = xx_project_info$archive
  )
  
  xx_about <- data.frame(
    bioproject = x[[1]]$Project$ProjectID$ArchiveID['accession'],
    id = x[[1]]$Project$ProjectID$ArchiveID['id'],
    name = xx_items$Name,
    title = xx_items$Title,
    descr = xx_items$Description
  )
  output <- list(
    about = xx_about,
    biosamples = xx_df,
    raw = xx_items
  )
  return(output)
}
eDetails <- list()
for(i in 1:length(eProject)){
  eDetails[[i]] <- extract_details(x = eProject[[i]])
}
eDetail_biosamples <- plyr::ldply(eDetails,function(x){x$biosamples})
bioproject_stats <- data.frame(table(eDetail_biosamples$bioproject))
bioprojects_stat <- bioproject_stats$Freq
names(bioprojects_stat) <- bioproject_stats$Var1
eDetail_about <- plyr::ldply(eDetails,function(x){x$about})
eDetail_about$nBiosamples <- bioprojects_stat[eDetail_about$bioproject]
eDetail_about <- eDetail_about[,c(1,ncol(eDetail_about),2:(ncol(eDetail_about)-1))]
write.csv(eDetail_about,
          file = 'processed-data/sra-00-bioprojects-statistics.csv',
          row.names = FALSE
)

write.csv(eDetail_biosamples,
          file = 'processed-data/sra-00-bioproject-with-biosamples.csv',
          row.names = FALSE
)
rm(list=ls())
# hello

