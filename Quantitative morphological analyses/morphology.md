Quantitative morphological analyses
================
BS
2022-12-08

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
data_raw    <- read.delim("data.txt", sep = "\t", header = T, check.names = F) 
conditions  <- read.delim("Conditions.txt", sep = "\t", header = T)
```

### log2 transform data

``` r
data_stat <- data_raw %>% 
      mutate(across(where(is.numeric), ~ log2(.x))) 
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

##### significance of clinical chemical parameter changes with 2 way anova

although p-values will be adjusted with two_way_anova_fn, only row
p-values will be used because these are separate/independent
measurements and multiple testing should not be an issue

``` r
anova_results <- two_way_anova_fn(data = data_stat, id_name = "Parameter", conditions_file = conditions, adjust_p_value = F, p_adj_method = "BH", add_l2fc = F)
```

    ## Joining, by = "Bioreplicate"
    ## `summarise()` has grouped output by 'Parameter'. You can override using the
    ## `.groups` argument.

### plot anova significant results (ratios)

#### function that obtains the data for error bar and plots bar diagrams

``` r
bar_chart_fn <- function(data_statistics, 
                         numeric_data,
                         conditions_data,
                         point_size,
                         id_name,
                         n_widht,
                         jitt_widht,
                         strip_text_size) {
  
  
  # prepare data for error bar calculation
  error_bar <- data_statistics %>% 
  left_join(numeric_data) %>% 
  select(contains(colnames(numeric_data))) %>% 
  pivot_longer(names_to = "Bioreplicate", values_to = "Value", -all_of(id_name)) %>% 
  left_join(conditions_data) %>% 
  group_by(!!as.symbol(id_name), Group) %>% 
  summarise(mean = mean(2^Value, na.rm=T), 
            sd = sd(2^Value, na.rm=T), 
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
  pivot_longer(names_to = "Bioreplicate", values_to = "Value", -all_of(id_name)) %>% 
  left_join(conditions_data) %>% 
  ungroup() %>% 
    left_join(error_bar) 
  
  # reorder data
  data_plot$Group <- factor(data_plot$Group, levels = c("NG", "HG"))
  
  
  # plot
  plot_data <- ggplot(data_plot %>% mutate(Value = 2^Value), aes(x=Group, y=Value))+
  geom_bar(stat = "summary", fun = mean, aes(fill = Group), alpha = 1, color = "black", lwd=0.5,
           width=n_widht/length(unique(data_plot$Group))) + 
  geom_jitter(size = point_size,  fill = "grey", alpha = 0.8, shape = 21, width = jitt_widht)+
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
  theme(strip.text = element_text(size = strip_text_size),
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
                                   filter(Parameter == "Relative liver weight"),
                              numeric_data = data_stat, point_size = 2.5,
                              n_widht = 0.7, jitt_widht = 0.05,
                              id_name = "Parameter", strip_text_size = 10, 
                              conditions_data = conditions)+
          ylab("")+
          geom_text(x = Inf, y = -Inf, 
          aes(label = paste("p=", `p-value group (HG/NG)`)),  
            size = 3.5, 
            hjust = 5, 
            #angle = 90,
            vjust = -25, 
            data = anova_results %>% 
            filter(Parameter == "Relative liver weight") %>% 
            mutate(across(where(is.numeric), round, 3)), inherit.aes = F)+
  facet_wrap(~Parameter, scales = "free_y", ncol = 3)+
  theme(legend.position = "bottom")+
  theme(plot.margin = margin(1,1,2,-3, "mm"))
```

    ## Joining, by = "Parameter"
    ## Joining, by = "Bioreplicate"
    ## `summarise()` has grouped output by 'Parameter'. You can override using the
    ## `.groups` argument.
    ## Joining, by = "Parameter"
    ## Joining, by = "Bioreplicate"
    ## Joining, by = c("Parameter", "Group")

![](morphology_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
ggsave("liver weight.svg", width = 3.5, height = 3.5)
```
