Statistical analysis
================
Sara Lindersson
2025-03-07

This script performs the statistical analysis and export the associated
figures for Lindersson et al 2025.

``` r
library(here)
library(tidyverse)
library(purrr)
library(ggplot2)
library(bbplot)
library(ggrepel)
library(broom)
library(ordbetareg)
library(tibble)
library(knitr)
library(forcats)
library(reshape2)
library(ggpubr)
library(ggradar)
library(fmsb)
```

## Load data

``` r
df <- readRDS("source_data.rds") %>%
  as.data.frame() %>%
  mutate(
    # Drop coordinates
    geometry = NULL)

variables <- c("uncertainty", "interdependency", "multi_levels", "volatility", "overlaps")
```

## How do scores vary across complexity dimensions?

``` r
# Create a long-format dataframe with variable names as a new column
df_long <- df %>%
  pivot_longer(cols = all_of(variables), 
               names_to = "dimension", 
               values_to = "value") %>%
  mutate(
    high_score = ifelse(value >= 4, T, F),
    low_score = ifelse(value <= 2, T, F)
  )

# Descriptive stats
kable(
  stats <-  df_long %>%
    group_by(dimension) %>%
    summarize(
      min = min(value),
      max = max(value),
      median = median(value),
      share_high_scores = round(mean(high_score), digits = 2),
      share_low_scores = round(mean(low_score), digits = 2),
      .groups = 'drop'
    ) %>%
    # Rearrange columns
    select(dimension, share_high_scores, share_low_scores, median, everything()) %>%
    arrange(desc(share_high_scores), desc(share_low_scores))
)

# Export stats
write.csv(stats, here("output", "tables", "1_stats_dimensions.csv"), row.names = FALSE)

# Calculate the proportion of each rating per dimension
df_proportions <- df_long %>%
  group_by(dimension, value) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(dimension) %>%
  mutate(proportion = count / sum(count)) 

hcl_palette <- hcl.colors(5, palette = "Emrld", rev = T)

fig <- df_proportions %>%
  ggplot(aes(x = proportion, y = fct_rev(factor(dimension)), fill = factor(value))) +
  geom_bar(stat = "identity") +
  scale_x_continuous(labels = scales::percent) +
  scale_fill_manual(values = hcl_palette, guide = guide_legend(title = NULL)) +
  labs(x = "Percentage of responses", y = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        legend.text = element_text(size = 8),
        axis.title = element_text(size = 8)
  ) +
  geom_text(aes(label = scales::percent(proportion, accuracy = 1)), 
            position = position_stack(vjust = 0.5),  # Position at the center of the bars
            size = 3, color = "white")  # Adjust size and color as needed

# Export figure
finalise_plot(plot_name = fig,
              save_filepath = here("output","figures","2_stacked_barplots_dimensions.pdf"),
              source_name = "",
              width_pixels = 350, height_pixels = 200)

# Do the ratings significantly vary across the dimensions?
# Kruskal-Wallis and Pairwise Wilcoxon
kable(
  kruskal_wilcoxon_score_dim <- bind_rows(
    tidy(kruskal.test(value ~ dimension, data = df_long)),
    tidy(pairwise.wilcox.test(df_long$value, df_long$dimension, p.adjust.method = "bonferroni")) %>%
      mutate(method = "Pairwise Wilcoxon test with Bonferroni")
  ) %>%
    # Round numeric columns
    mutate(
      statistic = round(statistic, digits = 2)
    ) %>%
    # Rearrange columns
    select(method, parameter, group1, group2, statistic, p.value)
)

# Export Kruskal-Wallis and Wilcoxon results
write.csv(kruskal_wilcoxon_score_dim, here("output", "tables", "2_kruskal_dimensions.csv"), row.names = FALSE)
```

## Co-occurrence of high- and low-scores across the dimensions

To identify if any dimensions seem to get high ratings together

``` r
# Function to compute Spearman correlation and p-values for all pairs
cor_and_pvalues <- function(df) {
  variables <- colnames(df)  # Get the names of the dimensions (columns)
  
  # Initialize empty matrices to store correlations and p-values
  cor_matrix <- matrix(NA, ncol = length(variables), nrow = length(variables), dimnames = list(variables, variables))
  pval_matrix <- matrix(NA, ncol = length(variables), nrow = length(variables), dimnames = list(variables, variables))
  
  # Compute Spearman correlation and p-value for each pair of columns
  for (i in 1:length(variables)) {
    for (j in i:length(variables)) {
      # Convert logical values to numeric (TRUE = 1, FALSE = 0)
      test <- cor.test(as.numeric(df[[variables[i]]]), as.numeric(df[[variables[j]]]), method = "spearman")
      
      cor_matrix[i, j] <- test$estimate  # Correlation coefficient
      pval_matrix[i, j] <- test$p.value  # P-value
      cor_matrix[j, i] <- test$estimate  # Mirror the result for symmetry
      pval_matrix[j, i] <- test$p.value  # Mirror the result for symmetry
    }
  }
  
  # Return correlation matrix and p-values
  list(correlation = cor_matrix, p_values = pval_matrix)
}

# Co-occurrence of high and low rankings
df_high <- df %>%
  mutate(across(all_of(variables), ~ .x >= 4)) %>%
  select(all_of(variables))

df_low <- df %>%
  mutate(across(all_of(variables), ~ .x <= 2)) %>%
  select(all_of(variables)) 

# Apply the function
result_high <- cor_and_pvalues(df_high)
result_low <- cor_and_pvalues(df_low)

# Create long dataframes for heatmaps
df_high_long <- melt(result_high$correlation) %>%
  rename(coef = value) %>%
  left_join(
    melt(result_high$p_values),
    by = c("Var1", "Var2")
  ) %>% 
  rename(p_value = value) %>%
  mutate(
    significance = factor(
      case_when(
        p_value <= 0.001 ~ 'P ≤ 0.001',
        p_value <= 0.01 ~ 'P ≤ 0.01, P > 0.001',
        p_value <= 0.05 ~ 'P ≤ 0.05, P > 0.01',
        T ~ 'P > 0.05'
      )
    )
  ) %>%
  mutate(
    Var1 = factor(Var1,
                  levels = c("interdependency", "multi_levels", "overlaps", "uncertainty", "volatility")),
    Var2 = factor(Var2,
                  levels = c("interdependency", "multi_levels", "overlaps", "uncertainty", "volatility"))
  )

df_low_long <- melt(result_low$correlation) %>%
  rename(coef = value) %>%
  left_join(
    melt(result_low$p_values),
    by = c("Var1", "Var2")
  ) %>% 
  rename(p_value = value) %>%
  mutate(
    significance = factor(
      case_when(
        p_value <= 0.001 ~ 'P ≤ 0.001',
        p_value <= 0.01 ~ 'P ≤ 0.01, P > 0.001',
        p_value <= 0.05 ~ 'P ≤ 0.05, P > 0.01',
        T ~ 'P > 0.05'
      )
    )
  ) %>%
  mutate(
    Var1 = factor(Var1,
                  levels = c("interdependency", "multi_levels", "overlaps", "uncertainty", "volatility")),
    Var2 = factor(Var2,
                  levels = c("interdependency", "multi_levels", "overlaps", "uncertainty", "volatility"))
  )

# Create heatmaps
fig_high <- ggplot(df_high_long, aes(Var1, fct_rev(Var2)), fill = "white") +
  geom_tile(color = "grey", fill = "white") +  # Add borders around tiles for clarity
  geom_point(aes(size = significance, color = abs(coef)), shape = 15, alpha = 0.8) +
  scale_color_gradientn(
    colors = hcl.colors(10, "Heat", rev = TRUE),
    limits = c(0, 0.5),
    name = NULL) +
  scale_size_manual(
    values = c(
      'P > 0.05' = 2,
      'P ≤ 0.05, P > 0.01' = 5,
      'P ≤ 0.01, P > 0.001' = 10,
      'P ≤ 0.001' = 15
    ),
    guide = guide_legend(),
    name = NULL
  ) +
  theme_minimal() +
  labs(title = "High rankings (4 and 5)", x = "", y = "") +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        legend.text = element_text(size = 8),
        axis.title = element_text(size = 8)
  ) +
  coord_fixed(ratio = 1)

fig_low <- ggplot(df_low_long, aes(Var1, fct_rev(Var2)), fill = "white") +
  geom_tile(color = "grey", fill = "white") +  # Add borders around tiles for clarity
  geom_point(aes(size = significance, color = abs(coef)), shape = 15, alpha = 0.8) +
  scale_color_gradientn(
    colors = hcl.colors(10, "Heat", rev = TRUE),
    limits = c(0, 0.5),
    name = NULL) +
  scale_size_manual(
    values = c(
      'P > 0.05' = 2,
      'P ≤ 0.05, P > 0.01' = 5,
      'P ≤ 0.01, P > 0.001' = 10,
      'P ≤ 0.001' = 15
    ),
    guide = guide_legend(),
    name = NULL
  ) +
  theme_minimal() +
  labs(title = "Low rankings (1 and 2)", x = "", y = "") +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        legend.text = element_text(size = 8),
        axis.title = element_text(size = 8)
  ) +
  coord_fixed(ratio = 1)

# Arrange figure
fig <- ggarrange(fig_low, fig_high, ncol = 2, nrow = 1, common.legend = T,
                 widths = c(1, 1), legend = "bottom")

# Export figure
finalise_plot(plot_name = fig,
              save_filepath = here("output","figures","3_heatmaps_dimensions.pdf"),
              source_name = "",
              width_pixels = 600, height_pixels = 400)
```

## How do scores vary between the hazard types?

``` r
# Remove all objects except the specified dataframes
rm(list = setdiff(ls(), c("df", "variables")))

# Create a long-format dataframe with hazard types as a new column
df_long <- df %>%
  pivot_longer(cols = all_of(variables), 
               names_to = "dimension", 
               values_to = "value") %>%
  mutate(
    high_score = ifelse(value >= 4, T, F),
    low_score = ifelse(value <= 2, T, F)
  )

# Descriptive stats
kable(
  stats <-  df_long %>%
    group_by(hazard) %>%
    summarize(
      n_cases = n()/5,
      min = min(value),
      max = max(value),
      median = median(value),
      share_high_scores = round(mean(high_score), digits = 2),
      share_low_scores = round(mean(low_score), digits = 2),
      .groups = 'drop'
    ) %>%
    # Rearrange columns
    select(hazard, n_cases, share_high_scores, share_low_scores, median, everything()) %>%
    arrange(desc(share_high_scores), desc(share_low_scores))
)

# Export stats
write.csv(stats, here("output", "tables", "3_stats_hazardtypes.csv"), row.names = FALSE)

# Calculate the proportion of each rating per hazard type
df_proportions <- df_long %>%
  group_by(hazard, value) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(hazard) %>%
  mutate(proportion = count / sum(count)) %>%
  mutate(hazard = factor(hazard, levels = c("earthquake", "volcanic activity", "drought", "multi-hazard", "flood", "epidemic", "storm", "wildfire")))

hcl_palette <- hcl.colors(5, palette = "Emrld", rev = T)

fig <- df_proportions %>%
  ggplot(aes(x = proportion, y = fct_rev(factor(hazard)), fill = factor(value))) +
  geom_bar(stat = "identity") +
  scale_x_continuous(labels = scales::percent) +
  scale_fill_manual(values = hcl_palette, guide = guide_legend(title = NULL)) +
  labs(x = "Percentage of responses", y = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        legend.text = element_text(size = 8),
        axis.title = element_text(size = 8)
  ) +
  geom_text(aes(label = scales::percent(proportion, accuracy = 1)), 
            position = position_stack(vjust = 0.5),  # Position at the center of the bars
            size = 3, color = "white")  # Adjust size and color as needed

# Export figure
finalise_plot(plot_name = fig,
              save_filepath = here("output","figures","4_stacked_barplots_hazards.pdf"),
              source_name = "",
              width_pixels = 350, height_pixels = 250)

# Kruskal-Wallis and Pairwise Wilcozon
kable(
  kruskal_wilcoxon_score_haz <- bind_rows(
    tidy(kruskal.test(value ~ hazard, data = df_long)),
    tidy(pairwise.wilcox.test(df_long$value, df_long$hazard, p.adjust.method = "bonferroni")) %>%
      mutate(method = "Pairwise Wilcoxon test with Bonferroni")
  ) %>%
    # Round numeric columns
    mutate(
      statistic = round(statistic, digits = 2)
    ) %>%
    # Rearrange columns
    select(method, parameter, group1, group2, statistic, p.value)
)

test <- kruskal.test(value ~ hazard, data = df_long)

# Export
write.csv(kruskal_wilcoxon_score_haz, here("output", "tables", "4_kruskal_hazards.csv"), row.names = FALSE)
```

## Do the complexity scores vary over time?

``` r
# Remove all objects except the specified dataframes
rm(list = setdiff(ls(), c("df", "variables")))

# Only include the cases with event years
df_y <- df %>% filter(!is.na(startyear_event)) 

### Spearman's rank

# Variables to test
variables_to_test <- c("uncertainty", "interdependency", "multi_levels",
                       "volatility", "overlaps", "score_median", "score_sum")

# Function to perform Spearman's rank correlation for a single variable
perform_spearman <- function(df, x_var, y_var) {
  test_result <- cor.test(
    x = df[[x_var]],
    y = df[[y_var]],
    method = "spearman"
  )
  
  # Tidy the result and add metadata
  tidy(test_result) %>%
    mutate(variable = x_var, tested_against = y_var)
}

# Apply the function to all variables and collect results
kable(
  spearman_results_time <- variables_to_test %>%
    map_dfr(~ perform_spearman(df_y, .x, "startyear_event")) %>%
    # Round numeric columns
    mutate(
      estimate = round(estimate, digits = 2),
      p.value = round(p.value, digits = 4)
    ) %>%
    # Rearrange columns
    select(method, alternative, variable, tested_against, everything()) %>%
    # Tidy variable names and order
    mutate(
      order = case_when(
        variable == "interdependency" ~ 1,
        variable == "multi_levels" ~ 2,
        variable == "overlaps" ~ 3,
        variable == "uncertainty" ~ 4,
        variable == "volatility" ~ 5,
        variable == "score_sum" ~ 7,
        variable == "score_median" ~ 6
      ),
      variable = case_when(
        variable == "interdependency" ~ "System interdependency",
        variable == "multi_levels" ~ "Multi-level governance",
        variable == "overlaps" ~ "Functional overlap",
        variable == "uncertainty" ~ "Uncertainty in data, models and predictions",
        variable == "volatility" ~ "Political volatility",
        variable == "score_sum" ~ "Total complexity score",
        variable == "score_median" ~ "Median complexity score"
      )
    ) %>% arrange(order) %>% select(-order)
)

# Export
write.csv(spearman_results_time, here("output", "tables", "5_spearman_years.csv"), row.names = FALSE)

manual_colors <- c("#DC8A5E","#FECC5B","#A8CADB","#16617B","#EA5559")

### Scatter plot of score_mean vs time
scatter <- ggplot(df_y, aes(y = score_sum, x = startyear_event)) +
  geom_smooth(method = "lm", se = TRUE, color = "darkgrey", fill = "lightgrey", lwd = 1) +
  geom_point(aes(color = hazard_group), size = 3.5, alpha = 0.7) +
  scale_color_manual(values = manual_colors) +
  labs(y = "Total complexity score across dimensions", x = "Start year of event") +
  scale_y_continuous(limits = c(0,30), breaks = seq(0, 30, by = 10)) +
  scale_x_continuous(limits = c(2003, 2025, breaks = seq(2005, 2025, by = 5))) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right") +
  guides(color = guide_legend(title = NULL))

# Export
finalise_plot(plot_name = scatter,
              save_filepath = here("output","figures","5_scatter_score_time.pdf"),
              source_name = "",
              width_pixels = 400, height_pixels = 250)
```

## How do the ratings vary across hazard types and dimensions?

``` r
# Remove all objects except the specified dataframes
rm(list = setdiff(ls(), c("df", "variables")))

variables <- c('interdependency', 'multi_levels', 'overlaps', 'uncertainty', 'volatility')

df2 <- df %>% select('case_id', 'hazard', all_of(variables))

earth <- ggradar(
  df2 %>% filter(hazard == 'earthquake') %>% select(-hazard),
  values.radar = c('1', '3', '5'),
  grid.min = 1, grid.mid = 3, grid.max = 5,
  group.line.width = .5,
  group.point.size = 1,
  fill = T,
  group.colours = "#FECC5B",
  background.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = 'none',
  axis.label.size = 4,
  grid.label.size = 4
) +
  labs(subtitle = 'Earthquake, N=4') +
  theme(plot.subtitle = element_text(size = 10)
  )

volcanic <- ggradar(
  df2 %>% filter(hazard == 'volcanic activity') %>% select(-hazard),
  values.radar = c('1', '3', '5'),
  grid.min = 1, grid.mid = 3, grid.max = 5,
  group.line.width = .5,
  group.point.size = 1,
  fill = T,
  group.colours = "#FECC5B",
  background.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = 'none',
  axis.label.size = 4,
  grid.label.size = 4
) +
  labs(subtitle = 'Volcanic activity, N=5') +
  theme(plot.subtitle = element_text(size = 10)
  )

drought <- ggradar(
  df2 %>% filter(hazard == 'drought') %>% select(-hazard),
  values.radar = c('1', '3', '5'),
  grid.min = 1, grid.mid = 3, grid.max = 5,
  group.line.width = .5,
  group.point.size = 1,
  fill = T,
  group.colours = "#DC8A5E",
  background.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = 'none',
  axis.label.size = 4,
  grid.label.size = 4
) +
  labs(subtitle = 'Drought, N=8') +
  theme(plot.subtitle = element_text(size = 10)
  )

multi <- ggradar(
  df2 %>% filter(hazard == 'multi-hazard') %>% select(-hazard),
  values.radar = c('1', '3', '5'),
  grid.min = 1, grid.mid = 3, grid.max = 5,
  group.line.width = .5,
  group.point.size = 1,
  fill = T,
  group.colours = "#EA5559",
  background.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = 'none',
  axis.label.size = 4,
  grid.label.size = 4
) +
  labs(subtitle = 'Multi-hazard, N=10') +
  theme(plot.subtitle = element_text(size = 10)
  )

flood <- ggradar(
  df2 %>% filter(hazard == 'flood') %>% select(-hazard),
  values.radar = c('1', '3', '5'),
  grid.min = 1, grid.mid = 3, grid.max = 5,
  group.line.width = .5,
  group.point.size = 1,
  fill = T,
  group.colours = "#16617B",
  background.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = 'none',
  axis.label.size = 4,
  grid.label.size = 4
) +
  labs(subtitle = 'Flood, N=13') +
  theme(plot.subtitle = element_text(size = 10)
  )

epidemic <- ggradar(
  df2 %>% filter(hazard == 'epidemic') %>% select(-hazard),
  values.radar = c('1', '3', '5'),
  grid.min = 1, grid.mid = 3, grid.max = 5,
  group.line.width = .5,
  group.point.size = 1,
  fill = T,
  group.colours = "#A8CADB",
  background.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = 'none',
  axis.label.size = 4,
  grid.label.size = 4
) +
  labs(subtitle = 'Epidemic, N=7') +
  theme(plot.subtitle = element_text(size = 10)
  )

storm <- ggradar(
  df2 %>% filter(hazard == 'storm') %>% select(-hazard),
  values.radar = c('1', '3', '5'),
  grid.min = 1, grid.mid = 3, grid.max = 5,
  group.line.width = .5,
  group.point.size = 1,
  fill = T,
  group.colours = "#16617B",
  background.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = 'none',
  axis.label.size = 4,
  grid.label.size = 4
) +
  labs(subtitle = 'Storm, N=7') +
  theme(plot.subtitle = element_text(size = 10)
  )

wildfire <- ggradar(
  df2 %>% filter(hazard == 'wildfire') %>% select(-hazard),
  values.radar = c('1', '3', '5'),
  grid.min = 1, grid.mid = 3, grid.max = 5,
  group.line.width = .5,
  group.point.size = 1,
  fill = T,
  group.colours = "#DC8A5E",
  background.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = 'none',
  axis.label.size = 4,
  grid.label.size = 4
) +
  labs(subtitle = 'Wildfire, N=4') +
  theme(plot.subtitle = element_text(size = 10)
  )

figure <- ggarrange(drought, wildfire, earth, volcanic,
                     epidemic, flood, storm, multi,
                     ncol = 4, nrow = 2)

# Export
finalise_plot(plot_name = figure,
              save_filepath = here("output","figures","6_radar_plots.pdf"),
              source_name = "",
              width_pixels = 600, height_pixels = 400)
```

End of script.
