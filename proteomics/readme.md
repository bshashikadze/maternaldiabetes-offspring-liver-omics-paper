statistical analysis of proteomics data and visualization
================
BS
23/11/2022

Dataset required for this script can be found at:
1. Liver_DIA_precursors_3d_old.tsv - can be downloaded from PRIDE repository (dataset ID - )

## load libraries

``` r
library(tidyverse)
library(randomForest)
library(missForest)
library(ggrepel)
library(ComplexHeatmap)
library(msigdbr)
library(circlize)
library(cowplot)
library(WebGestaltR)
library(ggpubr)
library(xlsx)
```

## Convert precursor data to protein data

### load data (DIA-NN main output)

``` r
precursors_diann_3d_old <- read.delim("Liver_DIA_precursors_3d_old.tsv", sep = "\t", header = T) 
```

### data filtering

filter for significance, signal quality, contaminants, unique peptides
adapted from:
<https://github.com/bshashikadze/diabetic_lung_proteomics_lipidomics/blob/main/DIA-NN%20quant%20with%20msempire/DIA-NN%20quant%20with%20msempire.Rmd>

``` r
diann_cleanup_function <- function(data, Q_Val, Global_Q_Val, Global_PG_Q_Val, Lib_Q_Val, Lib_PG_Q_Val, 
                                   experimental_library, unique_peptides_only, Quant.Qual, remove_contaminants) {
  
  # experimental library based analysis (e.g. GPF)
  if (experimental_library == T) {
  data <- data %>% 
  dplyr::filter(Q.Value           <= Q_Val) %>% 
  dplyr::filter(Global.Q.Value    <= Global_Q_Val) %>%
  dplyr::filter(Global.PG.Q.Value <= Global_PG_Q_Val)
  }
  
  # library free analysis (with mbr enabled)
  else {
    data <- data %>% 
    dplyr::filter(Q.Value         <= Q_Val)%>%
    dplyr::filter(Lib.Q.Value     <= Lib_Q_Val) %>% 
    dplyr::filter(Lib.PG.Q.Value  <= Lib_PG_Q_Val)
  }
  
  # filter for unique peptides and signal quality
  if (unique_peptides_only == TRUE) {
    unique <- 1
    }
  else {unique <- 0}
  
   data <- data %>%
   dplyr::filter(Proteotypic      >= unique) %>% 
   dplyr::filter(Quantity.Quality >= Quant.Qual) 
  
  # removing contaminant entries (according to maxquant common contaminants fasta file, which can be included during DIA-NN search)
   if (remove_contaminants == TRUE) {
  require(seqinr)
  contamintants           <- seqinr::read.fasta("contaminants.fasta")
  contaminant.names       <- seqinr::getName(contamintants) 
  data <- data %>%
  dplyr::filter(!stringr::str_detect(Protein.Group, stringr::str_c(contaminant.names, collapse="|")))
   return(data)
  }
  else {return(data)}
  }

# apply function on the data
data_filtered <- diann_cleanup_function(precursors_diann_3d_old, 
                                        Q_Val = 0.01, 
                                        Global_Q_Val = 0.01, 
                                        Global_PG_Q_Val = 0.01, 
                                        Lib_Q_Val = 0.01, 
                                        Lib_PG_Q_Val = 0.01, 
                                        experimental_library = T, 
                                        unique_peptides_only = TRUE, 
                                        Quant.Qual = 0.5, 
                                        remove_contaminants = T)
```

### from precursor to peptide, from peptide to protein

adapted from:
<https://github.com/bshashikadze/diabetic_lung_proteomics_lipidomics/blob/main/DIA-NN%20quant%20with%20msempire/DIA-NN%20quant%20with%20msempire.Rmd>

``` r
precursor_to_protein_function <- function(data, id_column, second_id_column, quantity_column, sum_charge){

  input   <- data
  results <- list()
  
  # subset for necessary columns
  data <- data %>% 
  dplyr::select("Precursor.Quantity", "Precursor.Normalised", "Run", all_of(quantity_column),
                "Stripped.Sequence", all_of(id_column)) 
  
  # aggregate charge states
  data_peptide <- data %>% 
  dplyr::select(-all_of(quantity_column)) 
        
    if (sum_charge == TRUE) {
  # summerise charge states by taking the sum
  data_peptide <- data_peptide %>% 
  dplyr::group_by(Run, Stripped.Sequence, !!as.symbol(id_column)) %>% 
  dplyr::summarise_all(sum, na.rm = TRUE) %>% 
  dplyr::ungroup()
      }
    else {
   
  # aggregate charge states by taking the precursor with the highest intensity  
  data_peptide <- data_peptide %>%
  dplyr::group_by(Run, Stripped.Sequence, !!as.symbol(id_column)) %>% 
  dplyr::summarise(Precursor.Quantity  = max(Precursor.Quantity), 
                   Precursor.Normalised = max(Precursor.Normalised)) %>% 
  dplyr::ungroup()
      }
 
  # from long to wide also import additional columns from the original data 
  data_peptide <- data_peptide %>% 
  tidyr::pivot_wider(names_from = "Run", 
                  values_from = c("Precursor.Quantity", "Precursor.Normalised"), 
                  c(all_of(id_column), "Stripped.Sequence"))

  data_peptide$peptide_q_value   <- input$Global.Q.Value[match(data_peptide[[id_column]], input[[id_column]])]
    
    
  # re-order
  data_peptide <- data_peptide %>% 
  dplyr::select(all_of(id_column), Stripped.Sequence,
                    starts_with("Precursor.Quantity"), starts_with("Precursor.Quantity"), 
                    starts_with("Precursor.Normalised"), peptide_q_value) 
    
  # remove entries without gene names
  if ("Genes" %in% colnames(data_peptide)) {
    data_peptide      <- data_peptide[data_peptide$Genes != "", ] 
  }
   
  # write results
  write.table(data_peptide, "peptides.txt", row.names = F, sep = "\t") 
  
  # protein level data
  data_pg <- data %>% 
  dplyr::select(Run, all_of(id_column), all_of(quantity_column)) %>%
  dplyr::distinct(Run, !!as.symbol(id_column), .keep_all = T) %>% 
  dplyr::mutate(Run = str_c("LFQ.intensity_", Run)) %>% 
  tidyr::pivot_wider(names_from = "Run", values_from = all_of(quantity_column), all_of(id_column))
  
  # count number of peptides
  n_pep <- data_peptide %>% 
    dplyr::group_by(!!as.symbol(id_column)) %>% 
    dplyr::summarise(n_pep = n_distinct(Stripped.Sequence)) 
  
  # match other columns
  n_pep$pg_Q_Val <- input$Global.PG.Q.Value[match(n_pep[[id_column]], input[[id_column]])]

  
  # get unique protein descriptions (distinct protein names for the same genes will be aggregated in one row separated by semicolon)
  protein_description <- input %>% 
  dplyr::select(!!as.symbol(id_column), First.Protein.Description) %>% 
  dplyr::distinct(!!as.symbol(id_column), First.Protein.Description, .keep_all = T) %>% 
  dplyr::group_by(!!as.symbol(id_column)) %>% 
  dplyr::summarize(First.Protein.Description=paste(First.Protein.Description,collapse=";")) %>% 
  dplyr::ungroup()


  # get unique protein descriptions (distinct protein names for the same genes will be aggregated in one row separated by semicolon)
  second_id <- input %>% 
  dplyr::select(!!as.symbol(id_column), !!as.symbol(second_id_column)) %>% 
  dplyr::distinct(!!as.symbol(id_column), !!as.symbol(second_id_column), .keep_all = T) %>%
  group_by(!!as.symbol(id_column)) %>% 
  dplyr::summarize(second_id =paste(!!as.symbol(second_id_column), collapse=";")) %>% 
  dplyr::ungroup()
  
  combined_additional <- protein_description %>% 
    dplyr::left_join(second_id)
  
  
  # remove entries without gene names
  if ("Genes" %in% colnames(data_peptide)) {
    data_pg      <- data_pg[data_pg$Genes != "", ] 
    }

  # join number of peptides and q-values
   data_pg <- data_pg %>% 
     dplyr::left_join(n_pep) %>% 
     left_join(combined_additional) 
  
  # write results
  write.table(data_pg, "proteingroups.txt", row.names = F, sep = "\t") 
  
  results[[1]] <- data_peptide
  results[[2]] <- data_pg

  return(results)}

# apply function on the data
peptide_and_pg_data <- precursor_to_protein_function(data_filtered, 
                                                 id_column = "Genes", 
                                                 second_id_column="Protein.Group", 
                                                 quantity_column = "Genes.MaxLFQ.Unique", 
                                                 sum_charge = TRUE)

rm(data_filtered, precursors_diann_3d_old, precursor_to_protein_function, diann_cleanup_function)
```

## load conditions file

``` r
conditions <- read.delim("Conditions.txt", sep = "\t", header = T)
```

## missing value imputation (random forest)

### data preparation

#### reoder columns (consistent with other data)

``` r
# protein groups data
proteingroups  <- peptide_and_pg_data[[2]] 

# reorder, first Genes, afterwards samples (from NG -> HG as in conditions file) finally additional columns
proteingroups_reordered <- proteingroups[c("Genes", str_c("LFQ.intensity_", conditions$Bioreplicate), colnames(proteingroups)[(nrow(conditions)+2):ncol(proteingroups)])]

# prepared protein groups for the suppl table and save
proteingroups_reordered_suppl <- proteingroups_reordered %>% 
  mutate(Intensity = rowSums(select(., starts_with("LFQ.")), na.rm = T)) %>% 
  arrange(desc(Intensity)) %>% 
  select(Genes, First.Protein.Description, second_id, starts_with("LFQ."), n_pep, pg_Q_Val) %>% 
  rename("Protein names"   = First.Protein.Description,
         "Protein groups"  = second_id,
         "Unique peptides" = n_pep,
         "q-value"         = pg_Q_Val) %>%
          write.table("proteingroups.txt", sep = "\t", row.names = F, quote = F)
```

#### count percentage of missing values

``` r
# data numeric
proteingroups  <- proteingroups_reordered %>% 
  column_to_rownames("Genes") %>% 
  select(starts_with("LFQ.")) %>% 
  mutate_all(log2)

# count missingness
n_total   <- nrow(proteingroups) * ncol(proteingroups)
n_missing <- sum(colSums(is.na(proteingroups)))

perc_missing <- (n_total/n_missing)
cat(paste0(round(perc_missing), " %"), "of the data is missing")
```

    ## 18 % of the data is missing

``` r
rm(n_total, n_missing, perc_missing)
```

#### filter for valid values (keep proteins with at least 60% valid values)

``` r
proteingroups_filtered <-  proteingroups %>% 
  mutate(n_valid = 100 - (100 * rowSums(is.na(.))/ncol(proteingroups))) %>% 
  filter(n_valid >= 60) %>% 
  select(-n_valid)

# count percentage of missing values after filtering
n_total_afterfiltering   <- nrow(proteingroups_filtered) * ncol(proteingroups_filtered)
n_missing_afterfiltering <- sum(colSums(is.na(proteingroups_filtered)))

perc_missing_afterfiltering <- (n_missing_afterfiltering/n_total_afterfiltering)
cat(paste0(perc_missing_afterfiltering, " %"), "of the data is missing after filtering")
```

    ## 0.0232595062316425 % of the data is missing after filtering

``` r
rm(n_total_afterfiltering, n_missing_afterfiltering, perc_missing_afterfiltering)
```

### missing value imputation

``` r
set.seed(12345)
system.time(data_imputed_RF <- missForest::missForest(as.matrix(proteingroups_filtered)))
```

    ##    user  system elapsed 
    ## 1511.91    2.72 1514.94

#### save imputed protein groups data

``` r
# clean-up the data
data_imputed <- data_imputed_RF[["ximp"]] %>% 
  as.data.frame() %>% 
  rownames_to_column("Genes") %>% 
  left_join(peptide_and_pg_data[[2]] %>% select(-starts_with("LFQ."))) 

# save results
write.table(data_imputed, "proteingroups_imputed.txt", sep = "\t", row.names = F)
```

## differential abundance analysis using 2 way ANOVA

### prepare data

``` r
# clean up proteomics data (1. retain only id column and columns containing protein intensities; 2. remove the string LFQ.intensity_ from each column, this is to be able to properly match with conditions file)
data_stat <- data_imputed %>% 
  select(Genes, starts_with("LFQ.")) %>% 
  rename_all(~str_replace_all(., 'LFQ.intensity_', ''))
```

### fold-change calculation

``` r
fc_function <- function(data, conditions_data, condition_name, compared_to, id_name, values_log) {
  # should values be log transformed?
  if (values_log == TRUE) {
    data_long <- data %>% 
    pivot_longer(names_to = "Bioreplicate", values_to = "Intensity", -!!as.symbol(id_name)) %>%
    left_join(conditions_data)
  }
  
  if (values_log == FALSE) {
    data_long <- data %>% 
    pivot_longer(names_to = "Bioreplicate", values_to = "Intensity", -!!as.symbol(id_name)) %>%
    mutate(Intensity = log2(Intensity)) %>% 
    left_join(conditions_data)
  }
  
  #fold-change calculation    
    data_fc <- data_long %>% 
    group_by(!!as.symbol(id_name), !!as.symbol(condition_name)) %>% 
    summarise(grp_mean = mean(Intensity, na.rm=T)) %>% 
    ungroup() 
  
    
  if (compared_to %in% unique(data_fc[[2]][2])) {
    data_fc <- data_fc %>% 
     select(-all_of(condition_name)) %>% 
     group_by(!!as.symbol(id_name)) %>% 
     mutate(l2fc=grp_mean-lag(grp_mean)) %>% 
     select(-grp_mean)  %>% 
     drop_na()
     }
  
   else{
    data_fc <- data_fc  %>% 
    select(-all_of(condition_name)) %>%  
    group_by(!!as.symbol(id_name))  %>% 
    mutate(l2fc=grp_mean-lead(grp_mean)) %>% 
    select(-grp_mean)   %>% 
    drop_na()
     }
    
    cat(paste("positive fold change means up in", compared_to, sep=" ")) 
    
    return(data_fc)
    
    }
# calculate fold-change separately for group and sex
l2fc_group <- fc_function(data_stat, condition = "Group", conditions_data = conditions,
                          compared_to  = "HG", values_log= T, id_name = "Genes")
```


``` r
l2fc_sex <- fc_function(data_stat, condition = "Sex", conditions_data = conditions,
                          compared_to  = "F", values_log= T, id_name = "Genes")
```

    ## Joining, by = "Bioreplicate"
    ## `summarise()` has grouped output by 'Genes'. You can override using the
    ## `.groups` argument.

    ## positive fold change means up in F

``` r
# merge fold-changes
l2fc_genes <- l2fc_group %>% 
  rename("l2fc group (HG/NG)" = l2fc) %>% 
  left_join(l2fc_sex) %>% 
  rename("l2fc sex (F/M)" = l2fc) 
```


### 2 way anova

#### perform anova analysis

##### define functions that performs 2 way anova

``` r
two_way_anova_fn <- function(data, id_name, conditions_file, adjust_p_value, p_adj_method, l2fc_data, add_l2fc) {
  
  data_anova  <- data %>% 
          pivot_longer(names_to = "Bioreplicate", 
          values_to = "Value", -all_of(id_name)) %>% 
          left_join(conditions_file) %>% 
          drop_na(Value) %>% 
          group_by(!!as.symbol(id_name)) %>% 
          summarise(`p-value` = 
          summary(aov(Value ~ Group*Sex))[[1]][["Pr(>F)"]][1:3]) 
  
  # correct all resulting p-values (pool) for multiple hypothesis testing
  if (adjust_p_value == TRUE) {
    
  data_anova$`Adjusted p-value`  <- p.adjust(data_anova$`p-value`, method = p_adj_method)
  # prepare empty data frame with proper comparisons
  anova_factors <- as.data.frame(rep(c("group (HG/NG)", "sex (F/M)", "group:sex"), length = nrow(data_anova)))
  # rename column 
  names(anova_factors) <- "Comparison"
  
  # final anova results
  anova_results  <- as.data.frame(cbind(data_anova, anova_factors)) %>% 
  pivot_wider(names_from = Comparison, 
                          values_from = c(`p-value`, `Adjusted p-value`), all_of(id_name), names_sep = " ") 
    
  }
  
  if (adjust_p_value == FALSE) {
  anova_factors <- as.data.frame(rep(c("p-value group (HG/NG)", "p-value sex (F/M)", 
                                       "p-value group:sex"), length = nrow(data_anova)))
  # rename column 
  names(anova_factors) <- "Comparison"
  
  # final anova results
  anova_results  <- as.data.frame(cbind(data_anova, anova_factors)) %>% 
  pivot_wider(names_from = Comparison, 
                          values_from = `p-value`, all_of(id_name), names_sep = " ") 
    
  }
  # add fold changes
  if (add_l2fc == TRUE) {
    anova_results  <- anova_results  %>% 
       left_join(l2fc_data)
  }
  return(anova_results)
  }
```

##### significance of protein changes with 2 way anova

``` r
anova_results <- two_way_anova_fn(data = data_stat, id_name = "Genes",                                  
                                  adjust_p_value = T,
                                  p_adj_method = "BH",
                                  add_l2fc = T,
                                  conditions_file = conditions, 
                                  l2fc=l2fc_genes)
```


#### clean-up anova results

``` r
# final anova results
anova_results  <- anova_results %>% 
  left_join(peptide_and_pg_data[[2]] %>% select(Genes, First.Protein.Description, second_id)) %>% 
  select("Genes", "First.Protein.Description", "second_id", "l2fc group (HG/NG)", "p-value group (HG/NG)", 
         "Adjusted p-value group (HG/NG)", "l2fc sex (F/M)", "p-value sex (F/M)", 
         "Adjusted p-value sex (F/M)", "p-value group:sex", "Adjusted p-value group:sex") %>% 
  arrange(-desc(`Adjusted p-value group (HG/NG)`)) %>% 
  rename(`Protein names` = First.Protein.Description,
         `Protein accession` = second_id) 


# separately significant factors and interactions

# significant by group
anova_results_group <- anova_results %>% 
  select(1:6) %>% 
  filter(`Adjusted p-value group (HG/NG)` <= 0.05 & abs(`l2fc group (HG/NG)`) > log2(1.5)) %>% 
          mutate(`Differentially abundant` = case_when(
                 `l2fc group (HG/NG)` > log2(1.5) ~  "increased in HG",
                 `l2fc group (HG/NG)` < -log2(1.5) ~ "decreased in HG",
             TRUE ~ "n.s."
           )) %>% 
  arrange(desc(`l2fc group (HG/NG)`)) 


# significant by sex
anova_results_sex <- anova_results %>% 
  select(1:3, 7:9) %>% 
  filter(`Adjusted p-value sex (F/M)` <= 0.05 & abs(`l2fc sex (F/M)`) > log2(1.5)) %>% 
          mutate(`Differentially abundant` = case_when(
                 `l2fc sex (F/M)` > log2(1.5) ~  "increased in F",
                 `l2fc sex (F/M)` < -log2(1.5) ~ "decreased in F",
             TRUE ~ "n.s."
           )) %>% 
  arrange(desc(`l2fc sex (F/M)`)) 


# interaction
anova_results_int <- anova_results %>% 
  select(1:3, 10,11) %>% 
  filter(`Adjusted p-value group:sex` <= 0.05) %>% 
  arrange(desc(`Adjusted p-value group:sex`))
```

#### define functions that performs 2 way anova Tukey’s honest significance difference (interaction significant proteins)

``` r
tkhsd_fn <- function(data, id_name, conditions_file, arrange_based, numeric_data) {
  
  # prepare data
  data_tukey <- data %>% 
              select(all_of(id_name)) %>%
              left_join(data_stat) %>% 
              pivot_longer(names_to = "Bioreplicate", values_to = "Intensity", 
                           -all_of(id_name)) %>% 
              left_join(conditions_file) %>% 
    drop_na(Intensity)
  
  # for loop. for each feature which showed significant interaction from
  # 2 way anova, THSD is calculated. significant pairs are extracted
  # make an empty list, where to each feature significant interactions will be assigned
 
   my_vec <- list()
  
  # for which feature
  
  sig_hits <- unique(data_tukey[[id_name]])
  
  # for loop
  for (i in sig_hits) {
    sub_df              <- data_tukey[data_tukey[[id_name]] %in% i,]
    anova_model         <- aov(data= sub_df,  Intensity ~ Group*Sex)
    anova_tukey         <- TukeyHSD(anova_model)
    tuk_interactions    <- anova_tukey[3][[1]] %>% 
    as.data.frame() %>% 
    rownames_to_column() %>% 
    filter(`p adj` < 0.05) %>% 
    select(rowname)
    tuk_interactions_t  <- t(tuk_interactions)
    sig_tuk_interaction <- matrix(apply(tuk_interactions_t,1,paste,collapse=";"),
                                nrow=1)
    my_vec[i] <- sig_tuk_interaction
  }
 
  # final data with anova and THSD statistics
  tkhsd <- list()
  tkhsd[[1]] <- anova_results_int %>% 
               left_join(my_vec %>% 
               unlist() %>% 
               as.data.frame() %>% 
               rownames_to_column("parameter") %>% 
               rename_all(~str_replace(., "parameter", id_name)) %>% 
               rename(`THSD pair` = 2)) %>% 
  arrange(desc(arrange_based))
  
  
  # data for interaction plot
  tkhsd[[2]] <- numeric_data %>%  
  left_join(tkhsd[[1]]) %>% 
  drop_na(`THSD pair`) %>% 
  select(1:ncol(data_stat)) %>% 
  pivot_longer(values_to = "Intensity", names_to = "Bioreplicate", -all_of(id_name)) %>% 
  left_join(conditions_file) %>% 
  mutate(joinedgr = str_c(Group, Sex, sep = "_")) %>% 
  group_by(!!as.symbol(id_name), joinedgr) %>% 
  summarise(mean = mean(Intensity, na.rm = T), 
            sd = sd(Intensity, na.rm = T), t.score = qt(p=0.05/2, 
            df=length(conditions_file$Bioreplicate),lower.tail=F), 
            ci = t.score * sd ) %>% 
  ungroup() %>% 
  separate(joinedgr, c("Group", "Sex")) 
  
  return(tkhsd)
}
```

##### THSD of protein interactions

``` r
anova_results_int_tuk <- tkhsd_fn(data = anova_results_int,  id_name = "Genes", 
                                  numeric_data = data_stat, 
                                  arrange_based = "Adjusted p-value group:sex", conditions_file = conditions)
```

    ## Joining, by = "Genes"
    ## Joining, by = "Bioreplicate"
    ## Joining, by = "Genes"
    ## Joining, by = "Genes"
    ## Joining, by = "Bioreplicate"
    ## `summarise()` has grouped output by 'Genes'. You can override using the
    ## `.groups` argument.

## volcano plot

### prepare data for the volcano plot

``` r
data_volcano <- anova_results %>% 
                mutate(significant = case_when(
                `Adjusted p-value group (HG/NG)` < 0.05 & `l2fc group (HG/NG)` > log2(1.5) | `Adjusted p-value group (HG/NG)` < 0.05 & 
                `l2fc group (HG/NG)` < -log2(1.5) ~ "+", TRUE ~ "n.s."),
                diff_abundant = case_when(
                `Adjusted p-value group (HG/NG)` < 0.05 & `l2fc group (HG/NG)` > log2(1.5) ~ "Increased_in_HG",
                `Adjusted p-value group (HG/NG)` < 0.05 & `l2fc group (HG/NG)` < -log2(1.5) ~ "Decreased_in_HG",
             TRUE ~ "n.s."
           )) 
```

### plot volcano plot

``` r
plot_volcano <- ggplot(data_volcano %>%                       
                       arrange(desc(diff_abundant)), 
                       mapping = aes(x = `l2fc group (HG/NG)`, y = -log10(`p-value group (HG/NG)`), 
                                     fill=diff_abundant, label = Genes, 
                                     alpha = diff_abundant))+
         geom_point(aes(shape =diff_abundant, size = diff_abundant), stroke = 0.25)+
         scale_shape_manual(values = c(n.s. = 16, Decreased_in_HG =21, Increased_in_HG =21))+
         scale_size_manual(values=c(n.s. = 1, Decreased_in_HG =1.6, Increased_in_HG =1.6))+
         scale_fill_manual(values=c("n.s." = "#999999", 
                                    "Decreased_in_HG"= "#0088AA", "Increased_in_HG"="#e95559ff"))+
         scale_alpha_manual(values= c("n.s." = 0.3, "Decreased_in_HG"= 1, "Increased_in_HG"= 1))+
         geom_text_repel(data = subset(data_volcano, 
                                       significant == "+" & `l2fc group (HG/NG)` > 0.6 | significant == "+"  & `l2fc group (HG/NG)` < -0.7),
                                       aes(label = Genes),
                                       size = 1.8,
                                       seed = 1234,
                                       color = "black",
                                       box.padding = 0.3,
                                       max.overlaps = 17,
                                       alpha = 1,
                          min.segment.length = 0)+
         theme_bw()+
         scale_x_continuous(limits = c(-2, 2.4)) +
         scale_y_continuous(limits = c(0,8), breaks = c(0,2.5,5,7.5)) +
         theme(panel.border = element_rect(size = 1, colour = "black"), 
               panel.grid.major = element_line(), 
               panel.grid.minor = element_blank(),
               panel.background = element_blank(), 
               axis.ticks = element_line(colour = "black"),
               axis.line = element_blank())+
                   theme(legend.position = "none", 
        legend.box.spacing = unit(0.8, 'mm'), 
        legend.title = element_blank(), 
        legend.text = element_blank())+
        theme(axis.title = element_text(size  = 9), 
               axis.text.x = element_text(size = 9, colour = "black", vjust = -0.1), 
               axis.text.y = element_text(size = 9, colour = "black"))+
         xlab("log2 fold change (HG/NG)")+
         ylab("-log10 p-value")
```

## plot anova interactions

``` r
# reorder data
anova_results_int_tuk[[2]]$Group <- factor(anova_results_int_tuk[[2]]$Group, levels = c("NG", "HG"))


# order of the facets (most significant first)
order <- anova_results_int %>% arrange(-desc(`Adjusted p-value group:sex`)) %>% 
                       select(Genes)
order <- order$Genes

# plot
ggplot(anova_results_int_tuk[[2]], aes(x=Group, y=mean)) +
  geom_line(size = 0.8, aes(group = Sex, color = Sex))+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                 position=position_dodge(0.05))+
  geom_point(aes(color = Sex))+
  scale_color_manual(values = c("#F0E442", "#E69F00"))+
  theme_bw()+
  theme(panel.border = element_rect(size = 1, colour = "black"),
       panel.grid.major = element_blank(), 
       panel.grid.minor = element_blank(), 
       panel.background = element_blank(), 
       axis.line = element_blank(),
  legend.position = "NONE") +
  theme(axis.title = element_text(size = 10), 
        axis.text.x = element_text(size= 9, colour = "black", vjust = -0.1), 
        axis.ticks = element_line(colour = "black"),
        axis.text.y = element_text(size = 9, colour = "black"))+
  theme(legend.position = c(0.6,0.05), legend.title = element_text(size = 8.5), axis.text = element_text(size = 8.5))+
  facet_wrap(~factor(Genes, levels=order), scales = "free_y", ncol = 5) +
  ylab("Protein intensity (log2)")+
  xlab("")+
  geom_text(x = Inf, y = -Inf, 
            aes(label = paste("p=",`Adjusted p-value group:sex`)),  
            size = 2.8, 
            angle = 90,
            hjust = -0.7, 
            vjust = -0.7, 
            data = anova_results_int_tuk[[1]] %>% 
                        mutate(across(where(is.numeric), round, 3)), inherit.aes = F)+
  theme(strip.background = element_blank(), strip.text = element_text())+
  guides(color=guide_legend(nrow=1, byrow=TRUE))
```

![](proteomics_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
ggsave("anova_interactions.svg", width = 7.1, height = 8.4)
```

## supervised clustering

### PCA

``` r
# pca function
pca_fn <- function(data, conditions_file, id_name) {
  
          # caclulate principal components of transposed dataframe
          PCA <- prcomp(t(data %>% column_to_rownames(id_name)))
          
          # export results for ggplot
          # according to https://www.youtube.com/watch?v=0Jp4gsfOLMs (StatQuest: PCA in R)
          pca.data <- data.frame(Bioreplicate=rownames(PCA$x),
          X=PCA$x[,1],
          Y=PCA$x[,2])
          pca.var <- PCA$sdev^2
          pca.var.per <- round(pca.var/sum(pca.var)*100,1)
          
          # export loadings
          loadings <- as.data.frame(PCA[["rotation"]]) %>% 
            select(1:2)
          
          # assign conditions
          pca.data <- pca.data %>% 
             left_join(conditions_file) 

          # save the data in a list
          pca_data      <- list()
          pca_data[[1]] <- pca.data
          pca_data[[2]] <- pca.var.per[1:2]
          pca_data[[3]] <- loadings
          
          return(pca_data)
          }

data_pca <- pca_fn(data_stat, conditions, id_name = "Genes")
```


#### PCA plot

``` r
plot_pca <- ggplot(data=data_pca[[1]], aes(x=X*-1, y=Y, fill= Group, label = ID))+
geom_point(size = 3, aes(shape = Sex), stroke = 0.25)+
scale_fill_manual(values= c('HG' = "#e95559ff", 'NG' = "#0088AA")) +
scale_shape_manual(values = c('F' = 21, 'M' = 22)) + 
xlab(paste("PC 1 - ", data_pca[[2]][1], "%", sep=""))+
ylab(paste("PC 2 - ", data_pca[[2]][2], "%", sep=""))+
geom_hline(yintercept = 0, linetype = "dashed")+
geom_vline(xintercept = 0, linetype = "dashed")+
theme_bw() + theme(panel.border = element_rect(size = 1, colour = "black"), 
                   axis.ticks = element_line(colour = "black"),
                   axis.title = element_text(size = 9, colour="black"), 
                   axis.text.x = element_text(size=9, colour="black", vjust = -0.1), 
                   axis.text.y = element_text(size = 9, colour="black"),
panel.grid.major = element_line(), panel.grid.minor = element_blank())+
theme(legend.title = element_text(colour="black", size=9))+
guides(shape = guide_legend(order = 2, override.aes = list(stroke = 1, shape  = c(0,1))),
       col = guide_legend(order = 1))+
  theme(legend.position = "top", 
        legend.box.spacing = unit(0.8, 'mm'), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 8.5))+
  guides(fill = guide_legend(override.aes=list(shape=21))) #https://github.com/tidyverse/ggplot2/issues/2322
```

### hierarchical clustering

#### data standardization

``` r
hm_prep_fn <- function(data, id_name) {
  
  # data standardization
  cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
  }
 
  data_HM <- as.data.frame(t(apply(data_stat %>% column_to_rownames(id_name), 1, cal_z_score)))

  # make sure that order of colnames matchs to the order of samples in conditions file
  order   <- conditions$Bioreplicate
  data_HM <- data_HM[order]

  # replace long names with short id names
  colnames(data_HM) <- conditions$ID
  
  return(data_HM)
  
}

data_HM <- hm_prep_fn(data = data_stat, id_name = "Genes")
rm(hm_prep_fn)
```

#### ploting heatmap

``` r
#makes a list with all necessary data for heatmap
hmap_all_data      <- list()
hmap_all_data[[1]] <- colorRamp2(c(-1.3, 0, 1.3), c("#0088AA",  "white", "#e95559ff"))
hmap_all_data[[2]] <- list('Group'  = c('NG' = "#0088AA", 'HG' = "#e95559ff"),
                           'Sex'     = c('F' = "#F0E442", 'M' = "#E69F00"))
hmap_all_data[[3]] <- HeatmapAnnotation(df = conditions %>% select(ID, Group, Sex) %>% column_to_rownames("ID"), 
                            show_legend    = F, 
                            which          = 'col',
                            annotation_name_gp   = gpar(fontsize = 8),
                            annotation_name_side = "left",
                            col                  = hmap_all_data[[2]],
                            height               = unit(0.7, "cm"),
                            simple_anno_size_adjust = TRUE)
# plotting
hmap <- Heatmap(as.matrix(data_HM),
                show_heatmap_legend        = F,
                row_dend_width             = unit(0.5, "cm"),
                column_dend_height         = unit(0.5, "cm"),
                show_column_names = T,
                show_row_names = F,
                column_names_side = "top",
                row_title = NULL,
                clustering_distance_columns = "euclidean",
                clustering_distance_rows    = "euclidean",
                clustering_method_rows      = "ward.D",
                clustering_method_columns   = "ward.D",
                use_raster = F,
                top_annotation              = hmap_all_data[[3]],
                width                       = unit(3.2, "in"),
                height                      = unit(4.72, "in"), 
                col                         = hmap_all_data[[1]], 
                border                      = TRUE,
                column_names_gp             = gpar(fontsize = 8)) 
ht <- draw(hmap)
```


``` r
#calculate actual plot size
w1 = ComplexHeatmap:::width(ht)
w1 = convertX(w1, "inch", valueOnly = TRUE)
h1 = ComplexHeatmap:::height(ht)
h1 = convertY(h1, "inch", valueOnly = TRUE)
c(w1, h1)
```


``` r
#for cowplot
plot_hmap = grid.grabExpr(draw(ht))
```

#### ploting legends separatelly

``` r
# color legend
hm_legend = grid.grabExpr(color_mapping_legend(hmap@matrix_color_mapping, plot = T, title = NULL, legend_direction = c("horizontal"),  title_gp = gpar(fontsize = 7, fontface = "bold"), param = list(at = c(-1,  1), labels = c("low", "high"), legend_width = unit(2, "cm"), labels_gp = gpar(fontsize=8)),  labels_gp = gpar(fontsize = 7)))

# legend for the group
Glycemia = grid.grabExpr(color_mapping_legend(hmap@top_annotation@anno_list[["Group"]]@color_mapping, nrow = 1, plot = T,  title_gp = gpar(fontsize = 7, fontface = "bold"),  labels_gp = gpar(fontsize = 7)))

# legend for the sex
Sex = grid.grabExpr(scale = 1, color_mapping_legend(hmap@top_annotation@anno_list[["Sex"]]@color_mapping,  nrow = 1, plot = T, title_gp = gpar(fontsize = 7, fontface = "bold"),  labels_gp = gpar(fontsize = 7)))

# combine legends
lg <- plot_grid(hm_legend, Glycemia, Sex, ncol = 3, rel_widths = c(1,1,1), rel_heights = c(1,1,1))
ggsave("legend.svg", width =2.8, height = 1)
```

## over-representation analysis (ORA)

#### ORA with WebGestalt

### perform ora

``` r
ora_for_cluster_function <- function(data, regulation, fc_threshold, fdr_threshold) {
  
  # regulation string to significance threshold
  if (regulation == "Decreased") {
    
    data_ora <- data %>% 
      filter(`l2fc group (HG/NG)` < -log2(fc_threshold) & `Adjusted p-value group (HG/NG)` <= fdr_threshold) %>% 
      select(Genes) %>% 
      as.list()
   
  }
  
  if (regulation == "Increased") {
       data_ora <- anova_results %>% 
      filter(`l2fc group (HG/NG)` > log2(fc_threshold) & `Adjusted p-value group (HG/NG)` <= fdr_threshold) %>% 
      select(Genes) %>% 
      as.list() 
  
       }  

  
  # perform ora
  set.seed(1234)
  outputDirectory <- getwd()
  enrichResult    <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
  enrichDatabase="geneontology_Biological_Process_noRedundant", interestGene=data_ora[[1]],
  interestGeneType="genesymbol", referenceSet = "genome",
  referenceGeneType="genesymbol", isOutput=TRUE,
  outputDirectory=outputDirectory, projectName=paste0(regulation), sigMethod = "fdr", fdrThr = 0.05, fdrMethod = "BH", minNum = 10, maxNum = 300)
  
  # read ora output
  ora_results <- read.delim(paste0("Project_", regulation, "/", "enrichment_results_", regulation, ".txt"))
  
  # tidy ora output
  ora_results <- ora_results %>% 
    mutate(Regulation = paste(regulation, "in HG"))
  
  return(ora_results)
}

# apply the function
increased_ora <- ora_for_cluster_function(anova_results, regulation = "Decreased", fc_threshold = 1.3, fdr_threshold = 0.05)
```

    ## Loading the functional categories...
    ## Loading the ID list...
    ## Loading the reference list...
    ## Summarizing the input ID list by GO Slim data...

    ## Performing the enrichment analysis...
    ## Begin affinity propagation...
    ## End affinity propagation...
    ## Begin weighted set cover...
    ## End weighted set cover...
    ## Generate the final report...
    ## Results can be found in the C:/Users/shashikadze/Documents/GitHub/maternaldiabetes-offspring-liver-omics-paper/proteomics/Project_Decreased!

``` r
decreased_ora <- ora_for_cluster_function(anova_results, regulation = "Increased", fc_threshold = 1.3, fdr_threshold = 0.05)
```

    ## Loading the functional categories...
    ## Loading the ID list...
    ## Loading the reference list...
    ## Summarizing the input ID list by GO Slim data...

    ## Performing the enrichment analysis...
    ## Begin affinity propagation...
    ## End affinity propagation...
    ## Begin weighted set cover...
    ## Remain is 0, ending weighted set cover
    ## Generate the final report...
    ## Results can be found in the C:/Users/shashikadze/Documents/GitHub/maternaldiabetes-offspring-liver-omics-paper/proteomics/Project_Increased!

``` r
# combine ora data
ora_data <- increased_ora %>% 
  bind_rows(decreased_ora) %>% 
  select(Regulation, geneSet, description, size, overlap, enrichmentRatio, pValue, FDR, userId) %>% 
  rename('GO id'                = geneSet, 
         'Biological process'   = description, 
         'Gene set size'        = size,
         'Gene count'           = overlap,
         'p-value'              = pValue,
         'Genes mapped'         = userId,
         'Fold enrichment'      = enrichmentRatio)

# save results
write.table(ora_data, "ora_results.txt", sep = "\t", row.names = F, quote = F)
```

### choose the processes that will be displayed on the plot

``` r
# replace zero fdr with lowest reported fdr (if any)
ora_data$FDR[ora_data$FDR == 0] <- min(ora_data$FDR[ora_data$FDR > 0])

# choose processes to be displayed on plot
ora_data_plot <- ora_data %>% 
  filter(Regulation == "Decreased in HG" &  `GO id` %in% c("GO:0006260", "GO:0032200", "GO:0006289", "GO:0055088", "GO:0006643")|
         Regulation == "Increased in HG" &  `GO id` %in% c("GO:0016051", "GO:0051188", "GO:0033865", "GO:0006575", "GO:0042180"))
```

## ORA dot plot

``` r
# facet labels
f_labels <- data.frame(Cluster = c("Decreased in HG", "Increased in HG"), label = c("Decreased in HG", "Increased in HG"))


# order rows based on enrichment score for SKM
ora_data_plot$`Biological process` <- factor(ora_data_plot$`Biological process`, levels = ora_data_plot$`Biological process`
                             [order(ora_data_plot$`Fold enrichment`)])
# plot
plot_ora <- ggplot(ora_data_plot, aes(x = `Fold enrichment`, y= `Biological process`, fill = -log10(FDR), size = `Gene count`)) +
 geom_point(shape = 21)+
 theme_bw() +
 theme(panel.border = element_rect(size = 1),
                            axis.text.x = element_text(colour = "black", size = 9),
                            axis.title  = element_text(size = 9),
                            axis.text.y = element_text(size = 9, colour = "black"))+
                      xlab("Enrichment score")+
                      ylab("")+
                      theme(plot.title = element_text(size = 8, hjust=0.5,
                                                      face = "bold")) +
                      theme(panel.border = element_rect(size = 1, colour = "black"),
                            panel.grid.major = element_line(),
                            axis.ticks = element_line(colour = "black"),
                            panel.grid.minor = element_blank())+
                      scale_size_continuous(name = "Count", range = c(1.8, 4.2))+
                      scale_fill_gradient(name = "-log10(FDR)", low = "#007ea7", high = "#003459")+
                      theme(legend.position = "bottom", 
                            legend.box.spacing = unit(0.8, 'mm'), 
                            legend.title = element_text(size = 8), 
                            legend.text = element_text(size = 8))+
                      theme(legend.key.size = unit(4, 'mm'))+
                      facet_grid(~Regulation, scales = "free")+
    geom_text(x = Inf, y = -Inf, aes(label = label),  size = 3, hjust = 1.05, vjust = -0.7, data = f_labels, inherit.aes = F)+
    theme(strip.background = element_blank(), strip.text = element_blank())
```

## save results

### combine plots

``` r
p0     <- rectGrob(width = 1, height = 1)
plot_1 <- ggarrange(plot_pca, plot_volcano, labels = c("B", "C"), font.label = list(size = 17),
          ncol = 1, nrow = 2, widths = c(7.1-w1,7.1-w1), heights = c(7.1-w1,7.1-w1))
```


``` r
plot_2 <- ggarrange(plot_hmap, NA, nrow = 2, heights = c(h1, 5.82-h1), widths = c(w1,w1))
```


``` r
plot_3 <- plot_grid(plot_2, plot_1, rel_widths = c(w1, 7.1-w1), rel_heights = c(h1+(5.82-h1), h1+(5.82-h1)), labels = c('A'), label_size = 17)
plot_4 <- plot_grid(plot_3, p0, ncol = 1, rel_widths = c(1,1), rel_heights = c(h1+(5.82-h1), 0.1))
plot_5 <- plot_grid(plot_4, plot_ora, ncol = 1, rel_widths = c(w1+7.1-w1, w1+7.1-w1), 
                    rel_heights = c(h1+0.5+(5.82-h1), 2.6), labels = c('','D'),  label_size = 17)
ggsave("proteomics_main.svg", width = 7.1, height = h1+0.5+(5.82-h1)+2.6)
```

### save data in a supplementary tables

``` r
if (!file.exists("Supplementary table 1.xlsx")) {
  
# anova results (all)
write.xlsx(as.data.frame(anova_results), file = "Supplementary table 1.xlsx", sheetName = "Suppl table 1B", 
  col.names = TRUE, row.names = FALSE, append = T)

# anova results (group)
write.xlsx(as.data.frame(anova_results_group), file = "Supplementary table 1.xlsx", sheetName = "Suppl table 1C", 
  col.names = TRUE, row.names = FALSE, append = T)

# anova results (ora)
write.xlsx(as.data.frame(ora_data), file = "Supplementary table 1.xlsx", sheetName = "Suppl table 1D", 
  col.names = TRUE, row.names = FALSE, append = T)

# anova results (sex)
write.xlsx(as.data.frame(anova_results_sex), file = "Supplementary table 1.xlsx", sheetName = "Suppl table 1E", 
  col.names = TRUE, row.names = FALSE, append = T)

# anova results (interaction) 
write.xlsx(as.data.frame(anova_results_int_tuk[[1]]), file = "Supplementary table 1.xlsx", sheetName = "Suppl table 1F", 
  col.names = TRUE, row.names = FALSE, append = T)

}
```

``` r
save.image()
```
