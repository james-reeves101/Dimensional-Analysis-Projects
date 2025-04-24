#Establish filenames
fileNames <- list.files(
  path = "Crick/maggie" ,
  pattern= " *.fcs",
  full.names = TRUE,
  ignore.case = TRUE
)
#Read in flowset
fset <- flowCore::read.flowSet(fileNames, truncate_max_range = FALSE)

colnames <- colnames(fset@frames[["export_D01 5131 CD3dep WLSM_Live.fcs"]]@exprs)
colnames <- colnames[7:20]
## identify optimum cofactor all channels
cofactors <-  flowVS::estParamFlowVS(fset, channels= colnames)
cofactorsdt <- data.table::data.table(colnames, cofactors)

exprn2 <- NULL
for (x in fileNames) {
  ff <- flowCore::read.FCS(x, emptyValue = FALSE, truncate_max_range = FALSE)
  this_data <- data.table::as.data.table(ff@exprs)
  names(this_data) <- as.vector(ff@parameters$desc)
  print(x)
  this_data %>%
    select(ends_with("Area"))->
    this_data 
  colnames(this_data)[  colnames(this_data) %in% c("[AF color 1] - Area",
                                                   "[AF color 2] - Area" )] = "[AF color] - Area" 
  this_data %>% 
    #add file name
    dplyr::mutate(file = basename(ff@description$FILENAME)) -> this_data
  exprn2 <- rbind(exprn2, this_data)
}
remove(this_data)


for (x in 3:16) {
  Spectre::do.asinh(exprn2, use.cols = x, cofactor = cofactors[[x - 2]]) -> exprn2
}

colnames(exprn2)[18:31] <- paste0(colnames(exprn2[,3:16]), "_asinh")
#FIX COL names
exprn2 <- colnames(exprn)[ colnames(exprn2)] 

exprn2 %>% 
  select(ends_with("_asinh"))->
  