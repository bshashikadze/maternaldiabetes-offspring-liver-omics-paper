weight and HOMA-IR
================
BS
10/05/2023

## load libraries

``` r
library(tidyverse)
library(ggpubr)
```

## load data

``` r
weight            <- read.delim("weight.txt", sep = "\t", header = T) 
conditions        <- read.delim("Conditions.txt", sep = "\t", header = T)
HOMA_IR           <- read.delim("HOMA_IR.txt", sep = "\t", header = T)
```

## differential abundance analysis using 2 way ANOVA

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
  anova_factors <- as.data.frame(rep(c("group (PHG/PNG)", "sex (F/M)", "group:sex"), length = nrow(data_anova)))
  # rename column 
  names(anova_factors) <- "Comparison"
  
  # final anova results
  anova_results  <- as.data.frame(cbind(data_anova, anova_factors)) %>% 
  pivot_wider(names_from = Comparison, 
                          values_from = c(`p-value`, `Adjusted p-value`), all_of(id_name), names_sep = " ") 
    
  }
  
  if (adjust_p_value == FALSE) {
  anova_factors <- as.data.frame(rep(c("p-value group (PHG/PNG)", "p-value sex (F/M)", 
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

##### significance of metabolite changes with 2 way anova

``` r
anova_results <- two_way_anova_fn(data = weight, id_name = "Parameter", conditions_file = conditions, adjust_p_value = F,  add_l2fc = F)
```

    ## Joining with `by = join_by(Bioreplicate)`

    ## Warning: Returning more (or less) than 1 row per `summarise()` group was deprecated in
    ## dplyr 1.1.0.
    ## ℹ Please use `reframe()` instead.
    ## ℹ When switching from `summarise()` to `reframe()`, remember that `reframe()`
    ##   always returns an ungrouped data frame and adjust accordingly.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## `summarise()` has grouped output by 'Parameter'. You can override using the
    ## `.groups` argument.

    ## Warning: Specifying the `id_cols` argument by position was deprecated in tidyr 1.3.0.
    ## ℹ Please explicitly name `id_cols`, like `id_cols = all_of(id_name)`.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

### plot results

#### define function which performs plotting

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
  summarise(mean = mean(Value, na.rm=T), 
            sd = sd(Value, na.rm=T), 
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
  data_plot$Group <- factor(data_plot$Group, levels = c("PNG", "PHG"))
  
  
  # plot
  plot_data <- ggplot(data_plot, aes(x=Group, y=Value))+
  geom_bar(stat = "summary", fun = mean, aes(fill = Group), alpha = 1, color = "black", lwd=0.15,
           width=n_widht/length(unique(data_plot$Group))) + 
  geom_jitter(size = point_size,  aes(shape = Sex), fill = "grey", alpha = 0.8, stroke =0.25, width = jitt_widht)+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05)) +
  scale_fill_manual(values=c('PNG' = "#0088AA", 'PHG' = "#e95559ff")) +
  scale_shape_manual(values = c('F' = 21, 'M' = 22)) + 
  theme_bw() +
  xlab("") +
  theme(panel.border = element_rect(linewidth  = 1, colour = "black"),
        axis.title = element_text(size = 9, colour = "black"),
        axis.text = element_text(size = 9, colour = "black"),
        axis.ticks = element_line(colour = "black"),
        legend.box.spacing = unit(0.8, 'mm'), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 8))+
  theme(strip.text = element_text(size = strip_text_size),
        strip.background = element_blank(),
        legend.margin=margin(0,0,0,0), 
        legend.box.margin=margin(-10,-10,-3,-10)) 
  
  
  return(plot_data)
 
 
}
```

#### plot data

``` r
# plot body weight
body_w <- bar_chart_fn(data_statistics = anova_results %>% filter(Parameter == "Body weight (kg)"),
                            numeric_data = weight, 
                            point_size = 2, 
                            strip_text_size = 9,
                            n_widht = 0.8, 
                            jitt_widht = 0.05,
                            id_name = "Parameter", 
                            conditions_data = conditions)+
            ylab("Body weight [kg]")+
  guides(shape = guide_legend(order = 2, override.aes = list(stroke = 1, shape  = c(1,0))),
         fill = guide_legend(order = 1))+
  theme(legend.position = "bottom")
```

    ## Joining with `by = join_by(Parameter)`
    ## Joining with `by = join_by(Bioreplicate)`
    ## `summarise()` has grouped output by 'Parameter'. You can override using the
    ## `.groups` argument.
    ## Joining with `by = join_by(Parameter)`
    ## Joining with `by = join_by(Bioreplicate)`
    ## Joining with `by = join_by(Parameter, Group)`

``` r
# plot liver weight
rel_liver_w <- bar_chart_fn(data_statistics = anova_results %>% filter(Parameter == "Rel_liver_weight"),
                            numeric_data = weight, 
                            point_size = 2, 
                            strip_text_size = 9,
                            n_widht = 0.8, 
                            jitt_widht = 0.05,
                            id_name = "Parameter", 
                            conditions_data = conditions)+
            ylab("Rel. liver weight")+
  guides(shape = guide_legend(order = 2, override.aes = list(stroke = 1, shape  = c(1,0))),
         fill = guide_legend(order = 1))+
  theme(legend.position = "bottom")
```

    ## Joining with `by = join_by(Parameter)`
    ## Joining with `by = join_by(Bioreplicate)`
    ## `summarise()` has grouped output by 'Parameter'. You can override using the
    ## `.groups` argument.
    ## Joining with `by = join_by(Parameter)`
    ## Joining with `by = join_by(Bioreplicate)`
    ## Joining with `by = join_by(Parameter, Group)`

``` r
# combine plots
ggarrange(body_w, rel_liver_w, common.legend = T, legend = "bottom", labels = c("A", "B"), font.label = list(size = 17))
```

![](general_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
ggsave("weights.svg", width = 4, height = 2.5)
```

## HOMA-IR (Homeostatic Model Assessment for Insulin Resistance)

``` r
# make grouping
HOMA_IR_tidy <- HOMA_IR %>% 
  mutate(Group_c = str_c(Group, " (", str_to_lower(Sex), ")"))


# error bars
error_bar_lipid <- HOMA_IR_tidy %>% 
  group_by(Group_c) %>% 
  summarise(mean = mean(HOMA_IR, na.rm=T), 
            sd = sd(HOMA_IR, na.rm=T)) %>% 
  ungroup()
  
# prepare data for plotting
  HOMA_IR_tidy <-  HOMA_IR_tidy %>% 
    left_join(error_bar_lipid) 
```

    ## Joining with `by = join_by(Group_c)`

``` r
# which comparisons will be displayed
IR_comparisons <- list(c("PNG (f)", "PHG (f)"), c("PNG (m)", "PHG (m)"),  c("PNG (f)", "PNG (m)"), c("PHG (f)", "PHG (m)"))


# make new group order 
HOMA_IR_tidy$Group_c <- factor(HOMA_IR_tidy$Group_c, levels = c("PNG (m)", "PHG (m)", "PNG (f)", "PHG (f)"))


# plot
ggplot(HOMA_IR_tidy, aes(x=Group_c, y=HOMA_IR, fill = Group_c))+
  geom_bar(stat = "summary", fun = mean, aes(fill = Group_c), alpha = 1, color = "black", lwd=0.15,
           width=2.1/length(unique(HOMA_IR_tidy$Group_c))) + 
  geom_jitter(size = 2,  aes(shape = Sex), fill = "grey", alpha = 0.8, stroke =0.25, width = 0.1)+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05)) +
  scale_shape_manual(values = c('F' = 21, 'M' = 22)) + 
  scale_fill_manual(values=c("PNG (f)" = "#0088AA", "PHG (f)" = "#e95559ff", "PNG (m)" = "#0088AA", "PHG (m)" = "#e95559ff")) +
  theme_bw() +
  ylab("") +
  xlab("") +
  ggtitle("HOMA-IR")+
  scale_y_continuous(breaks = c(0,1,2,3), limits = c(0,3))+
  theme(panel.border = element_rect(linewidth  = 1, colour = "black"),
        axis.title = element_text(size =10, colour = "black"),
        axis.text = element_text(size = 10, colour = "black"),
        axis.ticks = element_line(colour = "black"),
        plot.title = element_text(size = 10, hjust = 0.5),
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10)) +
   stat_compare_means(comparisons = IR_comparisons, method = "t.test", method.args = list(var.equal = TRUE),
                     size = 3, label = "p.format")+
                      theme(legend.position = "none")+
  theme(plot.margin = margin(1,1,-2,1, "mm"))
```

    ## [1] FALSE

![](general_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
  ggsave("homa_ir.pdf", width = 2.8, height = 2.45)
```
