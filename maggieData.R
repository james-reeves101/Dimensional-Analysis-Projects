#Establish filenames
fileNames <- list.files(
  path = "Crick/maggie" ,
  pattern= " *.fcs",
  full.names = TRUE,
  ignore.case = TRUE
)

exprn <- NULL
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
  exprn <- rbind(exprn, this_data)
}

remove(this_data)


#######Filter out CD3+ cells:
#exprn %>%
#  filter(`CD3 : FITC - Area_asinh` < )
#


to.asinh <- colnames(exprn)[3:16]
exprn <- Spectre::do.asinh(exprn, use.cols = to.asinh , cofactor = 200)

#Plot to 
ggplot(exprn, aes(x = `CD16 : AF700 - Area_asinh` , y = `CD56 : APC - Area_asinh` )) +
  geom_hex(bins=256)+
  scale_fill_continuous(type="viridis")+
  facet_wrap(~file)


#pivot to visualise asinh tranformation across markers asQC
exprn %>% 
  select(
    contains("_asinh"),file
  ) %>% 
  pivot_longer(cols = -file,
               names_to = "marker",
               values_to = "expression"
               ) |> ggplot(aes(x=expression))+
    geom_density()+
  xlim(-1,1.5)+
  facet_wrap(~marker) 


#write & zip csv files to save data.table post-transformation
fwrite(
  exprn,
  "Crick/maggie/posttransform.csv.gz" ,
  compress = "gzip"
)

#Select columns which have been transformed using asinh suffix
cluster_cols <- colnames(exprn)[grep("_asinh$", colnames(exprn))]

#Set seed for random values
set.seed(42)

#Run UMAP Algorithmn
exprn |>
  dplyr::select(all_of(cluster_cols)) |>
  uwot::umap(
    n_neighbors = 60, n_threads = (parallel::detectCores()-2),
    n_sgd_threads = (parallel::detectCores()-2), 
    verbose = TRUE
  ) -> umap_exprn
  
umap_exprn <- data.table::as.data.table(umap_exprn)
colnames(umap_exprn) <- c("UMAP_X", "UMAP_Y")

#cbind umap values to exprn data.table
exprn <- cbind(exprn,umap_exprn)


umap_exprn |> ggplot(aes(x=UMAP_X, y=UMAP_Y))+
  geom_hex(bins=256)+
  scale_fill_continuous(type="viridis")+
  coord_equal()


remove(umap_exprn)


####LABELLING MARKERS ON UMAP

marker <- "CD16 : AF700 - Area_asinh"

ColrMin <- quantile(exprn[[marker]], probs=c(0.01), na.rm = TRUE)

ColrMax <- quantile(exprn[[marker]], probs=c(0.995), na.rm = TRUE)
#Colour scheme
colour_scheme <- colorRampPalette(viridis::viridis_pal(option = "viridis")(50))



#ggplot function labelling by median expression per bin using summary function of bins + scale_fill 
exprn |> ggplot(aes(x=UMAP_X, y=UMAP_Y))+
  geom_hex(bins=256)+
  ggtitle(marker)+
  stat_summary_hex(aes(z= .data[[marker]]),
  fun= "median", bins = 256)+
  scale_fill_gradientn(
    colours = colour_scheme(50),
    limits = c(ColrMin, ColrMax), oob = scales::squish
  )+
  coord_equal()

p <- list()
for (marker in cluster_cols){
  ColrMin <- quantile(exprn[[marker]], probs=c(0.01), na.rm = TRUE)
  ColrMax <- quantile(exprn[[marker]], probs=c(0.995), na.rm = TRUE)
  #Colour scheme
  colour_scheme <- colorRampPalette(viridis::viridis_pal(option = "viridis")(50))
   exprn |> ggplot(aes(x=UMAP_X, y=UMAP_Y))+
    geom_hex(bins=256)+
    ggtitle(marker)+
    stat_summary_hex(aes(z= .data[[marker]]),
                     fun= "median", bins = 256)+
    scale_fill_gradientn(
      colours = colour_scheme(50),
      limits = c(ColrMin, ColrMax), oob = scales::squish
    )+
    coord_equal() -> p[[marker]]
} 

gridExtra::grid.arrange(grobs=p)

#tSNE:
exprn |>
  dplyr::select(all_of(cluster_cols)) |>
 Spectre::run.fitsne(
   dims= 2,
   theta = 0.5,
   intervals_per_integer = 1,
  ) -> tsne_exprn

#Broken
# exprn <- cbind(c(exprn, tsne_exprn$FItSNE_X, tsne_exprn$FItSNE_Y))


#plot tSNE without overlay
tsne_exprn |> ggplot(aes(x=FItSNE_X, y=FItSNE_Y))+
    geom_hex(bins=256)+
    scale_fill_continuous(type="viridis")+
    coord_equal()


#tSNE plot with marker overlay
tsne_exprn |> ggplot(aes(x=FItSNE_X, y=FItSNE_Y))+
  geom_hex(bins=256)+
  ggtitle(marker)+
  stat_summary_hex(aes(z= .data[[marker]]),
                   fun= "median", bins = 256)+
  scale_fill_gradientn(
    colours = colour_scheme(50),
    limits = c(ColrMin, ColrMax), oob = scales::squish
  )+
  coord_equal()





p <- list()
for (marker in cluster_cols){
  ColrMin <- quantile(exprn[[marker]], probs=c(0.01), na.rm = TRUE)
  ColrMax <- quantile(exprn[[marker]], probs=c(0.995), na.rm = TRUE)
  #Colour scheme
  colour_scheme <- colorRampPalette(viridis::viridis_pal(option = "viridis")(50))
  exprn |> ggplot(aes(x=FItSNE_X, y=FItSNE_Y))+
    geom_hex(bins=256)+
    ggtitle(marker)+
    stat_summary_hex(aes(z= .data[[marker]]),
                     fun= "median", bins = 256)+
    scale_fill_gradientn(
      colours = colour_scheme(50),
      limits = c(ColrMin, ColrMax), oob = scales::squish
    )+
    coord_equal() -> p[[marker]]
} 

gridExtra::grid.arrange(grobs=p)


###Phenograph 

d.set <- exprn %>% 
  select(all_of(cluster_cols)) %>% 
  as.matrix()

Rphenograph_out <- iCellR::Rphenograph(
  data = d.set,
  k = 15
)


exprn$Phenograph_cluster <- igraph::membership(Rphenograph_out[[2]])


 exprn |> ggplot(aes(UMAP_X, UMAP_Y, fill= as.factor(Phenograph_cluster)))+
   geom_hex(bins=256)+
   labs(fill="Phenograph clusters")
 
 
exprn |> ggplot(aes(UMAP_X, UMAP_Y, colour = .data[[marker]]))+
  geom_hex(bins=256)+
  stat_summary_hex(aes(z = .data[[marker]]),
                 fun = "median", bins = 256) +
  scale_fill_gradientn(colours = c(colour_scheme(50)),
    limits = c(ColrMin, ColrMax), oob = squish)


#Summarise expression per cluster in order to create heatmap 
exp <- Spectre::do.aggregate(
  exprn,
  use.cols = cluster_cols,
  by = "Phenograph_cluster"
)

#Heatmap plot - exports to WD
Spectre::make.pheatmap(
  dat=exp,
  sample.col = "Phenograph_cluster",
  plot.cols = cluster_cols,
  file.name = "pheatmap of markers_normalised.png",
  normalise = TRUE
)

#Load Example data to follow along:
exprn2 <- data.table::fread("Crick/data/flow-panel-data/Blood/mouse_blood_umap_Pheno.csv.gz")

table(exprn2$Phenograph_cluster) |>
  barplot(horiz=TRUE, las=1, xlab = "Count", ylab="Cluster")
  
exprn2 %>% 
  dplyr::filter(Phenograph_cluster %in% c(20,28)) %>% 
  ggplot(aes(CD4_asinh, CD8a_asinh)) +
  geom_hex(bins=256) +
  facet_wrap(file~Phenograph_cluster)


#Manually delete obsolete/ artifcat cluster
exprn2 %>% 
  dplyr::filter(Phenograph_cluster != 20) -> exprn2


#MANUALLY ASSIGN CLUSTERS
exprn2$Cluster_name <- NA

exprn2 |> mutate(Cluster_name = dplyr::case_when(
  Phenograph_cluster %in% c(2,11) ~ "CD4+ T Cells",
  Phenograph_cluster %in% c(8,9) ~ "CD8+ T Cells",
  TRUE ~ Cluster_name)) -> exprn2

exprn2 |> ggplot(aes(UMAP_X, UMAP_Y, fill = Cluster_name)) +
  geom_hex(bins=256)

exprn2 |> mutate(Cluster_name = dplyr::case_when(
  Phenograph_cluster %in% c(21:28) ~ "CD4+ CD8+ T Cells",
  Phenograph_cluster %in% c(13,16) ~ "Intermediate T Cells",
  TRUE ~ Cluster_name)) -> exprn2

exprn2 |> mutate(Cluster_name = dplyr::case_when(
  Phenograph_cluster %in% c(7) ~ "NK Cells",
  Phenograph_cluster %in% c(5) & NK11_asinh > 0.8 ~ "NK T Cells",
  Phenograph_cluster ==  5  & NK11_asinh <= 0.8 ~ "Unknown",
  TRUE ~ Cluster_name)) -> exprn2
  
exprn2 |> ggplot(aes(UMAP_X, UMAP_Y, fill = Cluster_name)) +
  geom_hex(bins=256)+
  labs(fill = "Cluster names") +
  scale_fill_manual(values = glasbey(n=12))


#Manual labelling of clusters via csv file then left_joining 
populations <- data.table::fread("Crick/data/flow-panel-data/Blood/Cell-populations-identified.csv")

exprn2 <- left_join(exprn2, populations, by= "Phenograph_cluster")

exprn2 |> ggplot(aes(UMAP_X, UMAP_Y, fill = Cell_pops)) +
  geom_hex(bins=256)+
  labs(fill = "Cluster names") +
  scale_fill_manual(values = glasbey(n=12))

####MEM 

MEM.input <- data.table::as.data.table(exprn2[,c(14:20)])

#add Phenograph clusters 
MEM.input$cluster <- exprn2$Phenograph_cluster

MEM.output <- cytoMEM::MEM(
  MEM.input,
  transform = FALSE,
  zero.ref = TRUE,
  scale.matrix = "arcsinh",
  scale.factor = 1
)


cytoMEM::build_heatmaps(
  MEM.output,
  cluster.MEM = "both",
  cluster.medians = "none",
  cluster.IQRs = "none",
  display.thresh = 2,
  output.files = TRUE,
  labels = TRUE,
  only.MEMheatmap = FALSE
)


