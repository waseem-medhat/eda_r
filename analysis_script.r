
# NOTE: The script uses the `## @knitr` syntax which marks the code as chunks
#       to be referred to in an Rmd file. This was made to prevent copying/
#       pasting code between 2 files and maximize reproducibility.


## @knitr packages

library(dplyr)
library(ggplot2)
library(gridExtra)
library(DataExplorer)
library(GGally)


## @knitr data

# read data
hep <- read.csv("indian_liver_patient.csv")

# clean data
hep <- hep %>%
  mutate(Diseased = ifelse(Dataset == 1, "Yes", "No")) %>%
  select(-Dataset) %>% 
  rename(Alanine_Aminotransferase = Alamine_Aminotransferase,
         Alkaline_Phosphatase     = Alkaline_Phosphotase,
         Total_Proteins           = Total_Protiens)


## @knitr datatypes

# structure
str(hep)


## @knitr missing

# missing counts and percent
number_missing <- sapply(hep, function(x) {sum(is.na(x))})
percent_missing <- sapply(hep, function(x) {round(mean(is.na(x)) * 100, 2)})

# wrap them in a data frame
missing <- data.frame(number_missing, percent_missing)
missing


## @knitr uni_histograms

# plot all numeric variables as histograms
plot_histogram(hep, nrow = 3, ncol = 3, geom_histogram_args = list(bins = 25))


## @knitr uni_ag_outliers

# get outlier values
outliers <- hep %>% 
  filter(Albumin_and_Globulin_Ratio > 2) %>% 
  pull(Albumin_and_Globulin_Ratio)

# table printing function
print_table <- function(x) {
  
  # takes a table object and prints it in a better form
  paste(names(x), ": ", x, "\n", sep = "")
  
}

# print outliers
cat(
  "Outliers\n",
  print_table(table(outliers))
)


## @knitr uni_log_histograms

# univariate histogram drawing function
hist_log <- function(variable_name) {
  
  ggplot(NULL, aes(hep[[variable_name]])) +
    geom_histogram() + 
    scale_x_log10() +
    labs(x = variable_name)
  
}

# log-transformed histograms
grid.arrange(hist_log("Alanine_Aminotransferase"),
             hist_log("Aspartate_Aminotransferase"),
             hist_log("Alkaline_Phosphatase"),
             hist_log("Direct_Bilirubin"),
             hist_log("Total_Bilirubin"),
             ncol = 2)


## @knitr uni_bili_tables

# most common bilirubin values
bdir_table <- sort(table(hep$Direct_Bilirubin), decreasing = TRUE)
btot_table <- sort(table(hep$Total_Bilirubin), decreasing = TRUE)


## @knitr uni_bdir_below_1

# get bilirubin levels < 1
bdir_below_1 <- hep %>% 
  filter(Direct_Bilirubin < 1) %>% 
  pull(Direct_Bilirubin)

# frequencies of bilirubin < 1
cat(
  "Direct Bilirubin Frequencies (< 1 mg/dL)\n",
  print_table(table(bdir_below_1))
)


## @knitr uni_histograms_more_bins

# histograms with more bins
plot_histogram(hep, nrow = 3, ncol = 3, geom_histogram_args = list(bins = 50))


## @knitr uni_albumin_table

# albumin frequencies
albumin_table <- sort(table(hep$Albumin), decreasing = TRUE)


# print most frequent albumin values
cat(
  "Albumin Frequencies (Top 10 Values)\n",
  print_table(albumin_table[1:10])
)


## @knitr uni_age_histogram

# age histogram by year
ggplot(hep, aes(Age)) +
  geom_histogram(binwidth = 1) +
  scale_x_continuous(breaks = seq(0, max(hep$Age), 5))


## @knitr uni_age_histogram_zoom

# subset ages for multiples of 5
ages_5 <- filter(hep, Age %% 5 == 0)

# zoom in on age
ggplot(hep, aes(Age)) +
  geom_histogram(binwidth = 1) +
  geom_histogram(data = ages_5, binwidth = 1, color = "white",
                 fill = "tomato2") + # colored layer for emphasis
  scale_x_continuous(breaks = seq(0, max(hep$Age), 5)) +
  coord_cartesian(xlim = c(38, 76))


## @knitr uni_bars

# gender proportions
ggplot(hep, aes(Gender, ..count../sum(..count..))) +
  geom_bar() +
  coord_flip()

# disease proportions
ggplot(hep, aes(Diseased, ..count../sum(..count..))) +
  geom_bar() +
  coord_flip()


## @knitr uni_freq

# print gender and disease frequencies
cat(
  "Gender\n",
  print_table(table(hep$Gender)),
  "\nDiseased\n",
  print_table(table(hep$Diseased))
)


## @knitr bi_correlations

# correlation matrix
plot_correlation(hep, type = "continuous",
                 cor_args = list(use = "pairwise.complete.obs"))


## @knitr ag_cor_alb

# albumin and a/g correlation
ag_cor_alb <- cor(
  hep$Albumin,
  hep$Albumin_and_Globulin_Ratio,
  use = "complete.obs"
)

## @knitr log_vars

# log-transform variables
hep_log <- hep %>% 
  mutate(log_alt  = log10(Alanine_Aminotransferase),
         log_ast  = log10(Aspartate_Aminotransferase),
         log_alp  = log10(Alkaline_Phosphatase),
         log_btot = log10(Total_Bilirubin),
         log_bdir = log10(Direct_Bilirubin))


## @knitr bi_scatterplots

# scatterplot matrix
hep_log %>% 
  select_if(is.numeric) %>% 
  select(-Alanine_Aminotransferase, -Aspartate_Aminotransferase,
         -Alkaline_Phosphatase, -Total_Bilirubin, -Direct_Bilirubin) %>% 
  ggpairs(
    lower = list(continuous = wrap("points", alpha = 0.1, stroke = 0))
  )

## @knitr ag_cor_bdir

ag_cor_bdir <- cor(
  hep_log$log_bdir,
  hep_log$Albumin_and_Globulin_Ratio,
  use = "complete.obs"
)

noout <- !(hep$Albumin_and_Globulin_Ratio %in% outliers) # "no outliers"

ag_cor_bdir_noout <- cor(
  hep_log$log_bdir[noout],
  hep_log$Albumin_and_Globulin_Ratio[noout],
  use = "complete.obs"
)

## @knitr bi_ag_bdir_scatterplot

# scatterplot of a/g ratio against log direct bilirubin
ggplot(hep_log, aes(log_bdir, Albumin_and_Globulin_Ratio)) +
  geom_jitter(width = 0.05) +
  geom_smooth(method = "lm", se = FALSE)


## @knitr bi_bars

# gender/disease stacked bars
ggplot(hep, aes(Diseased, fill = Gender)) +
  geom_bar(position = "fill") +
  coord_flip()


## @knitr crosstab

hep %>% 
  select(Diseased, Gender) %>% 
  table()


## @knitr bi_boxplots

plot_boxplot(hep, by = "Diseased", nrow = 3, ncol = 3)


## @knitr bi_log_boxplots

# boxplots by disease
hep_log %>% 
  select(starts_with("log"), Diseased) %>% 
  plot_boxplot(by = "Diseased", nrow = 3, ncol = 2)


## @knitr bi_means

# get means grouped by disease
means <- hep %>% 
  group_by(Diseased) %>% 
  select_if(is.numeric) %>% 
  summarize_all(.funs = function(x) {mean(x, na.rm = TRUE)}) %>% 
  select(-Diseased) %>% 
  t()

# add column names
colnames(means) <- c("No", "Yes")

# print as data frame
data.frame(means)


## @knitr bi_age_gender_facet

# age histogram faceted by gender
ggplot(hep, aes(Age)) +
  geom_histogram(binwidth = 1) +
  geom_histogram(data = ages_5, binwidth = 1, color = "white",
                 fill = "tomato2") + 
  scale_x_continuous(breaks = seq(0, max(hep$Age), 5)) +
  facet_grid(rows = vars(Gender), scales = "free_y")


## @knitr bi_age_gender_zoom

# zoom in again
ggplot(hep, aes(Age)) +
  geom_histogram(binwidth = 1) +
  geom_histogram(data = ages_5, binwidth = 1, color = "white",
                 fill = "tomato2") + 
  scale_x_continuous(breaks = seq(0, max(hep$Age), 5)) +
  facet_grid(rows = vars(Gender), scales = "free_y") +
  coord_cartesian(xlim = c(38, 76))


# # # Multivariate # # #


## @knitr multi_disease_gender

# multivariate boxplot drawing function
gbox_multivar <- function(variable_name, log_scale = FALSE) {
  
  # take variable name as a string and create a boxplot
  # grouped by disease, colored by gender
  
  plt <- ggplot(hep, aes(Diseased, hep[[variable_name]], fill = Gender)) +
    geom_boxplot() +
    coord_flip() +
    labs(y = variable_name)
  
  if (log_scale) {
    plt + scale_y_log10()
  } else {
    plt
  }
  
}

# boxplots
grid.arrange(gbox_multivar("Age"),
             gbox_multivar("Albumin"),
             gbox_multivar("Albumin_and_Globulin_Ratio"),
             gbox_multivar("Total_Proteins"),
             gbox_multivar("Alanine_Aminotransferase"  , log_scale = TRUE),
             gbox_multivar("Aspartate_Aminotransferase", log_scale = TRUE),
             gbox_multivar("Alkaline_Phosphatase"      , log_scale = TRUE),
             gbox_multivar("Direct_Bilirubin"          , log_scale = TRUE),
             gbox_multivar("Total_Bilirubin"           , log_scale = TRUE),
             ncol = 2)


## @knitr multi_bdir_female_no

hep %>% 
  filter(Gender == "Female", Diseased == "No") %>% 
  pull(Direct_Bilirubin) %>% 
  summary()


## @knitr multi_age

ggplot(hep, aes(Age)) + 
  geom_histogram(binwidth = 1) +
  geom_histogram(data = ages_5, binwidth = 1, color = "white",
                 fill = "tomato2") +
  facet_grid(Gender ~ Diseased, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, max(hep$Age), 5)) +
  coord_cartesian(xlim = c(38, 76))


## @knitr multi_outliers

# scatterplot of A/G ratio vs. direct bilirubin
# colored by disease category
ggplot(hep_log, aes(Albumin_and_Globulin_Ratio, Alkaline_Phosphatase)) +
  geom_jitter(width = 0.05, size = 3, alpha = 0.4, stroke = 0) +
  scale_color_brewer(palette = "Set1") +  # brighter color for clarity
  facet_grid(Diseased ~ Gender)


## @knitr final_1_summary

# function for alt summary stats
summarize_alt <- function(df) {
  
  summarize(
    df,
    mean  = mean(Alanine_Aminotransferase),
    sd    = sd(Alanine_Aminotransferase),
    n     = n(),
    se    = sd/sqrt(n),
    ci_lo = t.test(Alanine_Aminotransferase, mu = 0)$conf.int[1],
    ci_hi = t.test(Alanine_Aminotransferase, mu = 0)$conf.int[2]
  )

}

# alt summary stats by disease
alt_summary <- hep %>% 
  group_by(Diseased) %>% 
  summarize_alt()
  

## @knitr final_1

alt_summary %>% 
  
  # recreate disease category names
  mutate(Diseased = ifelse(Diseased == "Yes", "Diseased", "Healthy")) %>% 
  
  # bars
  ggplot(aes(Diseased, mean, fill = "")) +
  geom_col(width = 0.6) +
  
  # create confidence interval error bar
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.2, color = "grey50") +
  
  # make the plot horizontal
  coord_flip() +
  
  # remove grey background
  theme_minimal() +
  
  # set colors and hide legend
  scale_fill_viridis_d(option = "E", guide = FALSE) +
  
  # labels
  labs(
    y = "Alanine Aminotransferase Level (IU/L)",
    x = "",
    title = paste("High Blood Levels of Liver Enzymes",
                  "are Associated with Liver Disease"),
    subtitle = "Alanine Aminotransferase is Higher in Liver Disease Patients",
    caption = "Error bars represent 95% confidence intervals."
  )
  

## @knitr final_2_summary

# alt summary stats by disease and gender
alt_summary_gender <- hep %>% 
  group_by(Diseased, Gender) %>% 
  summarize_alt()


## @knitr final_2

alt_summary_gender %>% 
  
  # recreate disease category names
  ungroup() %>% 
  mutate(Diseased = ifelse(Diseased == "Yes", "Diseased", "Healthy")) %>% 
   
  # bars
  ggplot(aes(Diseased, mean, fill = Gender)) +
  geom_col(position = "dodge") +
  
  # create confidence interval error bar
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), color = "grey50",
                width = 0.2, position = position_dodge(1)) +
  
  # make the plot horizontal
  coord_flip() +
  
  # remove grey background
  theme_minimal() +
  
  # set colors
  scale_fill_viridis_d(option = "E", direction = -1) +
  
  # labels
  labs(
    y = "Alanine Aminotransferase Level (IU/L)",
    x = "",
    title = paste("Males Might Be Affected by Liver Disease",
                  "More Severely than Females"),
    subtitle = paste("Alanine Aminotransferase is Higher",
                     "in Diseased Males vs. Females"),
    caption = "Error bars represent 95% confidence intervals."
  )


## @knitr final_3

# histogram
  ggplot(hep, aes(Age, fill = (Age %% 5 == 0))) +
  geom_histogram(binwidth = 1) +
  
  # breaks
  scale_x_continuous(breaks = seq(0, max(hep$Age), 5)) +
  
  # facet
  facet_grid(rows = vars(Gender), scales = "free_y", switch = "y") +
  
  # set colors and hide legend
  scale_fill_viridis_d(option = "E", guide = FALSE, direction = -1) +

  # remove grey background
  theme_minimal() +
  
  # labels
  labs(
    y = "Frequency",
    title = "Do Old Men Report Their Ages Rounded to the Nearest Five?",
    subtitle = paste("Ages that are Multiples of 5 are Common in Men",
                     "(> 40 Years Old)")
  )
