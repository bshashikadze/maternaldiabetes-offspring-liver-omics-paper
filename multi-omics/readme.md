multi-omics analysis
================
BS
09/12/2022

Datasets required for this script can be found at:
1. data_1_log2 - is output of the [script - line 301](https://github.com/bshashikadze/maternaldiabetes-offspring-liver-omics-paper/blob/main/proteomics/proteomics.Rmd)
2. data_2_log2 - online supplementary tables, supplementary table 2 (sheet 1)
3. [Conditions file](https://github.com/bshashikadze/maternaldiabetes-offspring-liver-omics-paper/blob/main/clinical%20chemistry/Conditions.txt) (with information about samples and groups)

script of co-inertia analysis and visualization was adapted from
<https://github.com/bshashikadze/diabetes-lung-omics-paper/tree/main/multiomics%20coinertia>

CIA is based on: Meng, C., B. Kuster, A.C. Culhane, and A.M. Gholami, A
multivariate approach to the integration of multi-omics datasets. BMC
Bioinformatics, 2014. 15(1): p. 162.
<https://bioconductor.org/packages/release/bioc/html/omicade4.html>

## load libraries

``` r
library     (tidyverse) 
library       (ggrepel)
library      (omicade4) 
library(MetabolAnalyze)
library        (ggpubr)
```

## prepare omics datasets

``` r
dir.create("input")
```

    ## Warning in dir.create("input"): 'input' already exists

### load files and tidy

``` r
# proteomics data
data_1_log2 <- read.delim("input/proteingroups_imputed.txt", sep = "\t") %>% 
  rename_all(~str_replace(., "LFQ.intensity_", "")) %>% 
  column_to_rownames("Genes") %>% 
  select(starts_with("D"))

# metabolomics  
data_2_log2 <- read.delim("input/metabolites.txt", sep = "\t") %>% 
  select(-Class) %>% 
  column_to_rownames("Metabolite") %>% 
  log2

# check if the order of samples are the same in all datasets
all(names(data_1_log2) == names(data_2_log2))
```

    ## [1] TRUE

``` r
# conditions
conditions  <- read.delim("Conditions.txt", sep = "\t")
```

### which omics data will be analysed?

``` r
#should correspond to the order to data_1_log2 and data_2_log2; add if there are more datasets
omics_data <- c("Proteomics", "Metabolomics") 
```

### creates additions folders

``` r
dir.create("tmp")
dir.create("output")
for (i in omics_data) {
dir.create(paste0("tmp/", i))  
}
```

### centering and unit scaling

Alternatively “unit” scaling, or “none”

``` r
center_and_scale_function <- function(data, which_omics, scaling = "unit") {
                    data <- data %>% 
                    t() %>% 
                    scale(center = T, scale = F) %>% 
                    MetabolAnalyze::scaling(scaling) %>% 
                    t() %>% 
                    as.data.frame() %>% 
                    rownames_to_column() %>% 
                    rename(Variable = rowname) %>% 
                    write_tsv(paste("tmp/", which_omics, "/", "dataMatrix.tsv", sep="")) %>% 
                    mutate(omics = which_omics)
                    variableMetadata <- data %>% 
                    select("Variable") %>% 
                    write_tsv(paste("tmp/", which_omics, "/", "variableMetadata.tsv", sep=""))
                    return(data)}
# apply the function to the data
data_1_log2_scaled  <-  center_and_scale_function(data_1_log2,  which_omics = "Proteomics")
data_2_log2_scaled  <-  center_and_scale_function(data_2_log2,  which_omics = "Metabolomics")
```

### generation of sampleMetadata.tsv with conditions information

``` r
generate_Conditions_function <- function() {
  data_1_log2 %>% 
  colnames() %>% 
  as.data.frame() %>% 
  rename(sampleMetadata=1) %>% 
  mutate(Condition = sampleMetadata) %>% 
  mutate(id = sampleMetadata) %>% 
  write_tsv(paste("tmp/", "sampleMetadata.tsv", sep=""))
  cat("sampleMetadata file was generated in the tmp folder. modify only Conditions column (second column) 
      and copy it in all omics folders without renaming")}
generate_Conditions_function()
```


### join omics data frames

``` r
# join is performed by rows
joined_data <- rbind(data_1_log2_scaled, data_2_log2_scaled)
```

## Co-inertial analysis

### prepare dataset suitable for co-inertia analysis

modified from:
<https://du-bii.github.io/module-6-Integrative-Bioinformatics/current/session4/dubii_session4_ProMetIS_practical.html>

``` r
# function that reads all omics and associated datasets (sample and variable metadata)
function_test    <- function(which_omics, data) {
     obj_name    <- paste0(which_omics, ".eset")
     obj_name    <- phenomis::reading(file.path(paste0("tmp/", which_omics, "/")), 
                                      report.c = "none", output.c = "set")
                    MultiDataSet::add_eset(data,
                                           obj_name,
                                           dataset.type = 
                                           paste0("tmp/", which_omics, "/log2trans_scaled/"),
                                           GRanges = NA)}
# makes an empty list to which omics datasets will be assigned;
# for loop is used to iterate between omics datasets
multi.mset <- MultiDataSet::createMultiDataSet() 
for (i in omics_data) {
  multi.mset <- function_test(which_omics = i, multi.mset)
}
# defines conditions (does not metter from which omics)
pdata.df <- Biobase::pData(multi.mset[[paste0("tmp/", omics_data[1], "/log2trans_scaled/")]])
Condition.vc <- pdata.df[,"Condition"]
print(table(Condition.vc)) 
```

    ## Condition.vc
    ## HG NG 
    ##  9 10

``` r
cat("check if Conditions are correctly assigned in the printed table above")
```


### co-inertia calculation

``` r
set.seed(1234)
Coinertia.data <- MultiDataSet::as.list(multi.mset)
Coinertia.plot <- omicade4::mcia(Coinertia.data,)
#visualize
omicade4::plot.mcia(Coinertia.plot, axes = 1:2,
                    phenovec = Condition.vc,
                    sample.lab = T, gene.nlab = 7, df.color=1:2)
```


``` r
# extract RV
RV.obs <- Coinertia.plot[["mcoa"]][["RV"]][2]
cat("The RV coefficient of this coinertia analysis is:", RV.obs[1])
```

    ## The RV coefficient of this coinertia analysis is: 0.8002614

## permutation based significance calculation

### define function which performs permutation based p-value calculation

group labels will be randomly permuted and RV coefficient will be
calculated for the permuted data (function works very slowly)

``` r
# created folder with permuted data
create_folder_function <- function() {
       omics_data_rand <- replace(omics_data, 
                                  omics_data==omics_data[1], 
                                  paste0(omics_data[1], "_randomized"))
       dir.create(paste0("tmp/", omics_data_rand[1]))
       my_files <- list.files(paste0("tmp/", omics_data[1]))
       my_files <- my_files[-which(my_files=="sampleMetadata.tsv")]
       file.copy(from = paste0("tmp/", omics_data[1], "/", my_files),   
                   to = paste0("tmp/", omics_data_rand[1],"/", my_files))}
# apply the function
create_folder_function()
```

``` r
# permutation
cia_permutation_function <- function() {
             omics_data_rand <- replace(omics_data, 
                                  omics_data==omics_data[1], 
                                  paste0(omics_data[1], "_randomized"))
             labels.rand <- pdata.df %>% 
                            rownames_to_column("sampleMetadata")
             expl.rand   <- labels.rand[sample (1:nrow(labels.rand)),]
                            write.table(expl.rand, 
                            paste0("tmp/", omics_data_rand[1], "/sampleMetadata.tsv"),
                            sep = "\t", 
                            quote  = F, row.names = F)
  multi.mset.random <- MultiDataSet::createMultiDataSet() 
  for (i in omics_data_rand) {multi.mset.random <- function_test(multi.mset.random, 
                                                            which_omics = i)}
  Coinertia.data.rand <- MultiDataSet::as.list(multi.mset.random)
  Coinertia.plot.rand <- omicade4::mcia(Coinertia.data.rand)
  RV <- Coinertia.plot.rand[["mcoa"]][["RV"]][2]
  return(RV)}
```

### apply the function

``` r
# replicate permutation function (n is number of permutations)
set.seed(12345)
n_perm = 500
RV.rand <- replicate (n = n_perm, cia_permutation_function())
# store RV-s
RVs <- c (RV.rand, RV.obs)
hist (RVs, nclass = 100)  
abline (v = RV.obs, col = 'red')  # red line to indicate where in the histogram is the observed value
```

![](multi-omics_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
# calculate p-value
P_RV_permuted <- sum (RVs >= RV.obs)/(n_perm + 1)
cat("Permutation based p-value is", P_RV_permuted[1])
```

    ## Permutation based p-value is 0.001996008

## extract data for plotting

### extract loadings data

``` r
loadings_for_plotting_function <- function() {
data <- cbind(joined_data, 
              as.data.frame(Coinertia.plot[["mcoa"]][["axis"]][["Axis1"]]),  
              as.data.frame(Coinertia.plot[["mcoa"]][["axis"]][["Axis2"]])) %>% 
              select(Variable, omics, `Coinertia.plot[["mcoa"]][["axis"]][["Axis1"]]`,   
                     `Coinertia.plot[["mcoa"]][["axis"]][["Axis2"]]`) %>% 
              dplyr::rename(loading1 = `Coinertia.plot[["mcoa"]][["axis"]][["Axis1"]]`) %>% 
              dplyr::rename(loading2 = `Coinertia.plot[["mcoa"]][["axis"]][["Axis2"]]`) 
return(data)}
# apply the function
data_loadings <- loadings_for_plotting_function() %>% 
  mutate(label = case_when(loading1 >  2.22 ~ Variable,
                           loading1 < -1.87 ~ Variable))
```

### extract components

``` r
components_for_plotting_function <- function() {
        data  <- as.data.frame(Coinertia.plot[["mcoa"]][["Tl1"]]) %>% 
                 rownames_to_column() %>% 
                 mutate(omics = str_extract(rowname, omics_data)) %>% 
                 mutate(omics = zoo::na.locf(omics)) %>% 
                 mutate(id = str_extract(rowname, pdata.df$id)) %>% 
                 left_join(pdata.df)
                 return(data)}
# apply the function
data_components <- components_for_plotting_function()
```


## plotting

### plot loading plot

``` r
# plotting
loadingplot <- ggplot(data_loadings %>%                               
arrange(desc(omics)), mapping = aes(x=loading1, y=loading2, color = omics, alpha = omics, label =label))+
geom_point(aes(shape = omics), size = 1.3)+
scale_color_manual(values=c('Proteomics' = "#999999", 'Metabolomics' = "black"))+
scale_shape_manual(values = c('Proteomics' = 16, 'Metabolomics' = 15))+
scale_alpha_manual(values= c('Proteomics' = 0.5, 'Metabolomics'= 1))+
geom_label_repel(box.padding = unit(0.2, "lines"),
                 point.padding = unit(0.2, "lines"),
                 max.overlaps = Inf,
                 alpha = 1, color = "black", size =2.2)+
xlab("Component 1")+
ylab("Component 2")+
geom_hline(yintercept = 0, linetype = "dashed")+
geom_vline(xintercept = 0, linetype = "dashed")+
scale_y_continuous(breaks = c(-2,-1,0,1,2))+
theme_bw() + theme(panel.border = element_rect(size = 1),
                 axis.title = element_text(size  = 9, colour="black"), 
                 axis.text.x = element_text(size = 9, colour="black"), 
                 axis.text.y = element_text(size = 9, colour="black"),
panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
theme(plot.margin = margin(5,1,1,1, "mm")) +
theme(legend.position = "top", legend.box.spacing = unit(0.5, 'mm'), 
      legend.title = element_blank(), 
      legend.text = element_text(size = 8))
```

### plot co-inertia

``` r
coinertiaplot <- ggplot(data_components, aes(x=Axis1, y=Axis2, color = Condition))+
geom_point(aes(shape = omics), size = 2.5) +
geom_segment(aes(x = Axis1_Proteomics, y = Axis2_Proteomics, 
                 xend = Axis1_Metabolomics, yend = Axis2_Metabolomics, color = Condition), size =0.5, 
                 data = data_components %>% 
                 pivot_wider(-rowname, names_from = omics, values_from = c(Axis1, Axis2)))+
scale_shape_manual(values = c('Proteomics' = 16, 'Metabolomics' = 15))+
scale_color_manual(values = c('HG' = "#e95559ff", 'NG' = "#0088AA"))+
xlab(paste("Co-inertia axis 1 - ", (round(Coinertia.plot$mcoa$pseudoeig[1]*100)), "%", sep=""))+
ylab(paste("Co-inertia axis 2 - ", (round(Coinertia.plot$mcoa$pseudoeig[2]*100)), "%", sep=""))+
geom_hline(yintercept = 0, linetype = "dashed")+
geom_vline(xintercept = 0, linetype = "dashed")+
theme_bw() + theme(panel.border = element_rect(size = 1),
                   axis.title = element_text(size  = 9, colour="black"), 
                   axis.text.x = element_text(size = 9, colour="black"), 
                   axis.text.y = element_text(size = 9, colour="black"),
panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
theme(plot.margin = margin(5,1,1,1, "mm")) +
annotate("text", x = 2, y = 2, size  = 3.3, label = paste0("RV=", round(RV.obs,2), ",", " p=", round(P_RV_permuted, 3)))+
theme(legend.position = "top", 
      legend.box.spacing = unit(0.5, 'mm'), 
      legend.title = element_blank(), 
      legend.text = element_text(size = 8))
```

## merge plots

``` r
ggarrange(coinertiaplot, loadingplot, widths = c(3.55, 3.55),
          heights = c(3.55,3.55),
          labels = c("A", "B"), font.label = list(size =17),
          ncol = 2, nrow = 1,  legend = "bottom",
          common.legend = T)
```

``` r
ggsave("output/multiomicsfigure.svg", width = 7.1, height = 3.55)
```
