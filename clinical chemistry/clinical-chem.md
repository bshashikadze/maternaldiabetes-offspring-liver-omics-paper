statistical analysis of clinical chemical paramters data with
visualization
================
BS
01/12/2022

## load libraries

``` r
library(tidyverse)
library(cowplot)
library(ggpubr)
library(xlsx)
```

    ## Warning: package 'xlsx' was built under R version 4.2.2

## load data

``` r
clinical_chem_data    <- read.delim("clinicalparams.txt", sep = "\t", header = T) 
conditions            <- read.delim("Conditions.txt", sep = "\t", header = T)
```

## differential abundance analysis using 2 way ANOVA

### log2 transform data

``` r
data_stat <- clinical_chem_data %>% 
      select(-Explanation) %>% 
      mutate(across(where(is.numeric), ~ log2(.x)))
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
                          compared_to  = "HG", values_log = T, id_name = "Parameter")
```

    ## Joining, by = "Bioreplicate"
    ## `summarise()` has grouped output by 'Parameter'. You can override using the
    ## `.groups` argument.

    ## positive fold change means up in HG

``` r
l2fc_sex <- fc_function(data_stat, condition = "Sex", conditions_data = conditions,
                          compared_to  = "F", values_log  = T, id_name = "Parameter")
```

    ## Joining, by = "Bioreplicate"
    ## `summarise()` has grouped output by 'Parameter'. You can override using the
    ## `.groups` argument.

    ## positive fold change means up in F

``` r
# merge fold-changes
l2fc_parameters <- l2fc_group %>% 
  rename("l2fc group (HG/NG)" = l2fc) %>% 
  left_join(l2fc_sex) %>% 
  rename("l2fc sex (F/M)" = l2fc) 
```

    ## Joining, by = "Parameter"

### 2 way anova

#### perform anova analysis

##### define functions that performs 2 way anova

``` r
two_way_anova_fn <- function(data, id_name, conditions_file, l2fc) {
  
  data_anova  <- data %>% 
          pivot_longer(names_to = "Bioreplicate", 
          values_to = "Intensity", -all_of(id_name)) %>% 
          left_join(conditions_file) %>% 
          drop_na(Intensity) %>% 
          group_by(!!as.symbol(id_name)) %>% 
          summarise(`p-value` = 
          summary(aov(Intensity ~ Group*Sex))[[1]][["Pr(>F)"]][1:3]) 
  # correct all resulting p-values (pool) for multiple hypothesis testing
  data_anova$`Adjusted p-value`  <- p.adjust(data_anova$`p-value`, method = "BH")
  # prepare empty data frame with proper comparisons
  anova_factors <- as.data.frame(rep(c("group (HG/NG)", "sex (F/M)", "group:sex"), length = nrow(data_anova)))
  # rename column 
  names(anova_factors) <- "Comparison"
  
  # final anova results
  anova_results  <- as.data.frame(cbind(data_anova, anova_factors)) %>% 
  pivot_wider(names_from = Comparison, 
                          values_from = c(`p-value`, `Adjusted p-value`), all_of(id_name), names_sep = " ") 
  
  # add fold changes
  if (exists("l2fc")) {
    
    anova_results  <- anova_results  %>% 
       left_join(l2fc)
  }
  return(anova_results)
  
  }
```

##### significance of metabolite changes with 2 way anova

although p-values will be adjusted with two_way_anova_fn, only row
p-values will be used because these are separate/independent
measurements and multiple testing should not be an issue

``` r
anova_results <- two_way_anova_fn(data = data_stat, id_name = "Parameter", conditions_file = conditions, l2fc=l2fc_parameters)
```

    ## Joining, by = "Bioreplicate"
    ## `summarise()` has grouped output by 'Parameter'. You can override using the
    ## `.groups` argument.
    ## Joining, by = "Parameter"

##### clean-up anova results

``` r
# final anova results
anova_results  <- anova_results %>% 
  left_join(clinical_chem_data %>% select(Parameter, Explanation)) %>% 
  select("Parameter", "Explanation", "l2fc group (HG/NG)", "p-value group (HG/NG)", 
         "l2fc sex (F/M)", "p-value sex (F/M)", 
         "p-value group:sex") %>% 
  arrange(-desc(`p-value group (HG/NG)`)) 


# separately significant factors and interactions
# significant by group
anova_results_group <- anova_results %>% 
  select(1:4) %>% 
  filter(`p-value group (HG/NG)` <= 0.05 & abs(`l2fc group (HG/NG)`) > log2(1)) %>% 
          mutate(`Differentially abundant` = case_when(
                 `l2fc group (HG/NG)` > log2(1) ~  "increased in HG",
                 `l2fc group (HG/NG)` < -log2(1) ~ "decreased in HG",
             TRUE ~ "n.s."
           )) %>% 
  arrange(desc(`l2fc group (HG/NG)`)) 


# significant by sex
anova_results_sex <- anova_results %>% 
  select(1:2, 5:6) %>% 
  filter(`p-value sex (F/M)` <= 0.05 & abs(`l2fc sex (F/M)`) > log2(1)) %>% 
          mutate(`Differentially abundant` = case_when(
                 `l2fc sex (F/M)` >  log2(1)  ~  "increased in F",
                 `l2fc sex (F/M)` < -log2(1)  ~  "decreased in F",
             TRUE ~ "n.s."
           )) %>% 
  arrange(desc(`l2fc sex (F/M)`)) 


# interaction
anova_results_int <- anova_results %>% 
  select(1:2, 7) %>% 
  filter(`p-value group:sex` <= 0.05) %>% 
  arrange(desc(`p-value group:sex`))
```

#### define functions that performs 2 way anova Tukey’s honest significance difference (interaction significant proteins)

``` r
tkhsd_fn <- function(data, id_name, conditions_file, filter_based, numeric_data) {
  
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
  arrange(filter_based)
  
  
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

##### THSD of metabolite interactions

``` r
anova_results_int_tuk <- tkhsd_fn(data = anova_results_int,  id_name = "Parameter", numeric_data = data_stat, filter_based = "p-value group:sex", conditions_file = conditions)
```

    ## Joining, by = "Parameter"
    ## Joining, by = "Bioreplicate"
    ## Joining, by = "Parameter"
    ## Joining, by = "Parameter"
    ## Joining, by = "Bioreplicate"
    ## `summarise()` has grouped output by 'Parameter'. You can override using the
    ## `.groups` argument.

## plot anova interactions

``` r
# reorder data
anova_results_int_tuk[[2]]$Group <- factor(anova_results_int_tuk[[2]]$Group, levels = c("NG", "HG"))

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
  theme(legend.position = "bottom", legend.title = element_text(size = 8.5), axis.text = element_text(size = 8.5))+
  facet_wrap(~factor(Parameter, levels=unique(anova_results_int_tuk[[1]]$Parameter)), scales = "free", ncol = 5) +
  ylab("Concentration (log2)")+
  geom_text(x = Inf, y = -Inf, 
            aes(label = paste("p=",`p-value group:sex`)),  
            size = 2.8, 
            angle = 90,
            hjust = -0.7, 
            vjust = -0.7, 
            data = anova_results_int_tuk[[1]] %>% 
                        mutate(across(where(is.numeric), round, 3)), inherit.aes = F)+
  theme(strip.background = element_blank(), strip.text = element_text())
```

![](clinical-chem_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
ggsave("anova_interactions_clinpar.svg", width = 3.5, height = 3.5)
```

### plot anova significant results (ratios)

#### function that obtains the data for error bar and plots bar diagrams

``` r
bar_chart_fn <- function(data_statistics, 
                         numeric_data,
                         conditions_data,
                         id_name,
                         n_widht,
                         jitt_widht) {
  
  
  # prepare data for error bar calculation
  error_bar <- data_statistics %>% 
  left_join(numeric_data) %>% 
  select(contains(colnames(numeric_data))) %>% 
  pivot_longer(names_to = "Bioreplicate", values_to = "Intensity", -all_of(id_name)) %>% 
  left_join(conditions_data) %>% 
  group_by(!!as.symbol(id_name), Group) %>% 
  summarise(mean = mean(2^Intensity, na.rm=T), 
            sd = sd(2^Intensity, na.rm=T), 
            n = n(),
            sem = sd/sqrt(n),
            t.score = qt(p=0.05/2, 
            df=length(conditions$Bioreplicate),lower.tail=F), 
            ci = t.score * sd ) %>% 
  ungroup()
  # prepare data for plotting
  data_plot <- data_statistics %>% 
  left_join(numeric_data) %>% 
  select(contains(colnames(numeric_data))) %>%  
  pivot_longer(names_to = "Bioreplicate", values_to = "Intensity", -all_of(id_name)) %>% 
  left_join(conditions_data) %>% 
  ungroup() %>% 
    left_join(error_bar) 
  
  # reorder data
  data_plot$Group <- factor(data_plot$Group, levels = c("NG", "HG"))
  
  
  # plot
  plot_data <- ggplot(data_plot %>% mutate(Intensity = 2^Intensity), aes(x=Group, y=Intensity))+
  geom_bar(stat = "summary", fun = mean, aes(fill = Group), alpha = 1, color = "black", lwd=0.5,
           width=n_widht/length(unique(data_plot$Group))) + 
  geom_jitter(size = 1.5,  fill = "grey", alpha = 0.8, shape = 21, width = jitt_widht)+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05)) +
  scale_fill_manual(values=c('NG' = "#0088AA", 'HG' = "#e95559ff")) +
  theme_bw() +
  xlab("") +
  theme(panel.border = element_rect(size = 1, colour = "black"),
        axis.title = element_text(size = 9, colour = "black"),
        axis.text = element_text(size = 9, colour = "black"),
        axis.ticks = element_line(colour = "black"),
        legend.box.spacing = unit(0.8, 'mm'), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 8))+
  theme(strip.text.y = element_text(size = 9), 
        strip.background = element_blank(),
        legend.margin=margin(0,0,0,0), 
        legend.box.margin=margin(-12,-12,-3,-12)) 
  
  
  return(plot_data)
 
}
```

#### plot data

``` r
# plot ratios
bar_chart_fn(data_statistics = anova_results %>% 
                                   filter(`p-value group (HG/NG)`<= 0.1),
                              numeric_data = data_stat,
                              n_widht = 0.7, jitt_widht = 0.05,
                              id_name = "Parameter", 
                              conditions_data = conditions)+
          ylab("")+
          geom_text(x = Inf, y = -Inf, 
          aes(label = paste("p=", `p-value group (HG/NG)`)),  
            size = 3, 
            hjust = -1.5, 
            angle = 90,
            vjust = -13.5, 
            data = anova_results %>% 
            filter(`p-value group (HG/NG)`<= 0.1) %>% 
            mutate(across(where(is.numeric), round, 3)), inherit.aes = F)+
  facet_wrap(~factor(Parameter, levels=c("Bilirubin (mg/dl)", "Albumin (g/dl)", 
                                         "NEFA (mmol/l)", "Glycerol (mmol/l)", "TG (mg/dl)", "ALT (U/l)")), scales = "free_y", ncol = 2)+
  theme(legend.position = "bottom")
```

    ## Joining, by = "Parameter"
    ## Joining, by = "Bioreplicate"
    ## `summarise()` has grouped output by 'Parameter'. You can override using the
    ## `.groups` argument.
    ## Joining, by = "Parameter"
    ## Joining, by = "Bioreplicate"
    ## Joining, by = c("Parameter", "Group")

    ## Warning: Removed 1 rows containing non-finite values (stat_summary).

    ## Warning: Removed 1 rows containing missing values (geom_point).

![](clinical-chem_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
ggsave("clinicalparameters.svg", width = 3.5, height = 5)
```

    ## Warning: Removed 1 rows containing non-finite values (stat_summary).
    ## Removed 1 rows containing missing values (geom_point).

### save data in a supplementary tables

``` r
if (!file.exists("Supplementary table 3.xlsx")) {
  
# metabolomics data (all)
write.xlsx(as.data.frame(clinical_chem_data), file = "Supplementary table 3.xlsx", sheetName = "Suppl table 1A", 
  col.names = TRUE, row.names = FALSE, append = T)
  
# anova results (all)
write.xlsx(as.data.frame(anova_results), file = "Supplementary table 3.xlsx", sheetName = "Suppl table 1B", 
  col.names = TRUE, row.names = FALSE, append = T)

# anova results (group)
write.xlsx(as.data.frame(anova_results_group), file = "Supplementary table 3.xlsx", sheetName = "Suppl table 1C", 
  col.names = TRUE, row.names = FALSE, append = T)

# anova results (sex)
write.xlsx(as.data.frame(anova_results_sex), file = "Supplementary table 3.xlsx", sheetName = "Suppl table 1D", 
  col.names = TRUE, row.names = FALSE, append = T)

# anova results (interaction) 
write.xlsx(as.data.frame(anova_results_int_tuk[[1]]), file = "Supplementary table 3.xlsx", sheetName = "Suppl table 1E", 
  col.names = TRUE, row.names = FALSE, append = T)
}
```

``` r
save.image()
```