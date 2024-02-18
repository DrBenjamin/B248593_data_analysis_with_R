## Setup 
install.packages("pacman")
library(pacman)
pacman::p_load(tidyverse, ggplot2, epitools, finalfit, webshot, broom, gt, here)

# Open Assessment instructions with Codebook section
browseURL(here("meta_data", "Assignment_instructions_2024.pdf"))


## Data input
framingham <- read_csv(here("raw_data", "framingham_2024.csv"))


## Data checking
View(framingham)
dim(framingham)
str(framingham)

# Showing total columns count
ncol(framingham)

# Showing total rows count
nrow(framingham)

# Showing columns with missing values
missing <- colSums(is.na(framingham))

# Just showing columns with missing values
col_missing <- missing[missing > 0]

# Writing column names with value in brackets, like cigsPerDay (32)
col_missing_names <- paste0(names(col_missing), " (", col_missing, ")")


## Data preparation
framingham_df <- framingham %>%
  mutate(
    sex = factor(sex, labels = c("male", "female")),
    currentSmoker = factor(currentSmoker, labels = c("No", "Yes")),
    BPMeds = factor(BPMeds, labels = c("No", "Yes")),
    prevalentStroke = factor(prevalentStroke, labels = c("No", "Yes")),
    prevalentHyp = factor(prevalentHyp, labels = c("No", "Yes")),
    diabetes = factor(diabetes, labels = c("No", "Yes")),
    anyCHD = factor(anyCHD, labels = c("No", "Yes"))
  )

# Viewing the dataframe
View(framingham_df)

# Saving the dataframe
write_csv(framingham_df, here("processed_data", "framingham_df.csv"))

# Cleaning up the environment
rm(framingham, missing, col_missing, col_missing_names)


## 1. We want to study the association between the occurrence of any coronary
## heart disease and two risk factors:
## a. What is the relative risk of occurrence of any coronary heart disease
## given the smoking status at the first examination?
# Contingency table
framingham_df %>%
  select(anyCHD, currentSmoker) %>%
  table()

# Checking assumption
subset_framingham <- framingham_df %>%
  drop_na(currentSmoker) %>%
  sample_n(1000)

# Distribution of currentSmoker
p0 <- subset_framingham %>%
  ggplot(aes(x = currentSmoker)) +
  geom_bar(aes(y = ..prop.., fill = anyCHD)) +
  labs(title = "Figure A1: Current Smoker",
       subtitle = "Distribution",
       x = "Current Smoker", y = "Proportion",
       caption = "Data source: Framingham Heart Study")
p0

# Saving the plot
ggsave(p0, file = here("figures/", "Plot_A1.png"))

# Relative risk
framingham_df %>%
  select(currentSmoker, anyCHD) %>%
  table() %>%
  epitab(., method = "riskratio", rev = "columns")

rm(subset_framingham)


## b. How is occurrence of any coronary heart disease associated with age category?
## Hint: The American Heart Association (AHA) reports that the incidence of
## cardiovascular disease (CVD) in US men and women is different between age
## groups 40–59 years, 60–79 years, and above the age of 80.
# Contingency table
framingham_df %>%
  select(anyCHD, age) %>%
  table()

# Range of age
framingham_df %>%
  select(age) %>%
  range()

# Mean and Median of the age variable
framingham_df %>%
  select(age) %>%
  summarise(
    age.mean = mean(age),
    age.median = median(age)
  )

# Creating crosstab
crosstab_age <- framingham_df %>%
  select(age, anyCHD)

# Contingency table
crosstab_age %>%
  table()

# Viewing crosstab
View(crosstab_age)

# Creating sample dataset
subset_age <-
  crosstab_age %>%
  drop_na(age) %>%
  sample_n(1000)

# Checking distribution of age (also for normality)
p1 <- subset_age %>%
  ggplot(aes(x = age)) +
  geom_histogram(aes(y = after_stat(density))) +
  geom_density(aes(y = after_stat(density))) + 
  labs(title = "Figure A2: Age",
       subtitle = "Distribution",
       x = "Age", y = "Density",
       caption = "Data source: Framingham Heart Study")
p1

# Saving the plot
ggsave(p1, file = here("figures/", "Plot_A2.png"))

# Shapiro-Wilk test
subset_age %>%
  pull(age) %>%
  shapiro.test()

# Showing histogram with of age in relation to `anyCHD`
p2 <- subset_age %>%
  ggplot(aes(x = age)) +
  geom_histogram(aes(y = after_stat(count), fill = ..density..), bins = 20) +
  labs(title = "Figure A3: Age distribution",
       subtitle = "in relation to the occurrence of any coronary heart disease",
       x = "Age", y = "Count",
       caption = "Data source: Framingham Heart Study") +
  facet_grid(rows = vars(anyCHD))
p2

# Saving the plot
ggsave(p2, file = here("figures/", "Plot_A3.png"))

# Summarising of age in relation to anyCHD
crosstab_age %>%
  summary_factorlist(
    dependent = "anyCHD", explanatory = "age",
    column = FALSE, total_col = TRUE
  )

# Performing t-test (is this parametric test valid?
# as age is NOT normally distributed!!!)
crosstab_age %>%
  t.test(age ~ anyCHD, data = .)

# Performing wilcox-test (non-parametric test, as age is NOT normally distributed)
# age in relation to anyCHD
crosstab_age %>%
  wilcox.test(age ~ anyCHD, data = ., conf.int = TRUE)

# Cleaning up the environment
rm(crosstab_age)

# Recoding age into 2 level factor (40-59, 60-79)
crosstab_age_cat <- framingham_df %>%
  filter(age >= 40) %>%
  mutate(age_cat = cut(age,
    breaks = c(39, 59, 79),
    labels = c("40-59", "60-79")
  )) %>%
  select(anyCHD, age_cat)

# Contingency table
crosstab_age_cat %>%
  table()

# Viewing crosstab
View(crosstab_age_cat)

# Checking age_cat
p3 <- crosstab_age_cat %>%
  mutate(anyCHD = fct_rev(anyCHD)) %>%
  ggplot(aes(x = age_cat, fill = anyCHD)) +
  geom_bar(position = "fill") +
  labs(title = "Figure A4: Age categories",
       subtitle = "in relation to the occurrence of any coronary heart disease",
       x = "Age category", y = "Proportion",
       caption = "Data source: Framingham Heart Study") +
  theme_minimal()
p3

# Saving the plot
ggsave(p3, file = here("figures/", "Plot_A4.png"))

# Relative risk
crosstab_age_cat %>%
  mutate(age_cat = fct_relevel(age_cat, "60-79")) %>%
  table() %>%
  epitab(., method = "riskratio", rev = "columns")

# Odds ratios
crosstab_age_cat %>%
  mutate(age_cat = fct_relevel(age_cat, "60-79")) %>%
  table() %>%
  epitab(., method = "oddsratio", rev = "columns")

# Cleaning up the environment
rm(crosstab_age, crosstab_age_cat, subset_age)


## 2. We want to study the association between risk factors:
## a. Is systolic blood pressure associated with smoking status?
# Contingency table
framingham_df %>%
  select(sysBP, currentSmoker) %>%
  table()

# Range of sysBP
framingham_df %>%
  select(sysBP) %>%
  range()

# Mean and Median of the sysBP variable
framingham_df %>%
  select(sysBP) %>%
  summarise(
    sysBP.mean = mean(sysBP),
    sysBP.median = median(sysBP)
  )

# Creating crosstab
crosstab <- framingham_df %>%
  select(currentSmoker, sysBP)

# Viewing crosstab
View(crosstab)

# Creating sampple dataset
subset_sysBP <- crosstab %>%
  drop_na(sysBP) %>%
  sample_n(1000)

# Checking distribution of sysBP (also for normality)
p4 <- subset_sysBP %>%
  ggplot(aes(x = sysBP)) +
  geom_histogram(aes(y = after_stat(density))) +
  geom_density(aes(y = after_stat(density))) +
  labs(title = "Figure A5: systolic blood pressure",
       subtitle = "Distribution in relation to smoking status",
       x = "Systolic blood pressure (mm Hg)", y = "Denstity",
       caption = "Data source: Framingham Heart Study") +
  facet_grid(rows = vars(currentSmoker))
p4

# Saving the plot
ggsave(p4, file = here("figures/", "Plot_A5.png"))

# Performing Shapiro-Wilk test
subset_sysBP %>%
  pull(sysBP) %>%
  shapiro.test()

# Summarising of smoking status in relation to sysBP
crosstab %>%
  summary_factorlist(
    dependent = "sysBP", explanatory = "currentSmoker",
    column = FALSE, total_col = TRUE
  )

# Performing Chi-squared test
crosstab %>%
  table() %>%
  chisq.test()

# Performing Wilcoxon rank sum test (non-parametric test,
# as `sysBP` is NOT normally distributed)
crosstab %>%
  wilcox.test(sysBP ~ currentSmoker, data = ., conf.int = TRUE)

# Lineare Regression
crosstab %>%
  lm(sysBP ~ currentSmoker, data = .) %>%
  summary()

# Cleaning up the environment
rm(crosstab)


## b. Is prevalence of hypertension associated with smoking status?
# Contingency table
framingham_df %>%
  select(prevalentHyp, currentSmoker) %>%
  table()

# Creating crosstab
crosstab <- framingham_df %>%
  select(currentSmoker, prevalentHyp)

# Viewing crosstab
View(crosstab)

# Summarising of smoking status in relation to sysBP
crosstab %>%
  summary_factorlist(
    dependent = "prevalentHyp", explanatory = "currentSmoker",
    column = FALSE, total_col = TRUE
  )

# Checking distribution of prevalentHyp
p5 <- crosstab %>%
  mutate(prevalentHyp = fct_rev(prevalentHyp)) %>%
  ggplot(aes(x = currentSmoker, fill = prevalentHyp)) +
  geom_bar(position = "fill") +
  labs(title = "Figure A6: Prevalence of hypertension",
       subtitle = "in relation to smoking status",
       x = "Smoking status", y = "Proportion",
       caption = "Data source: Framingham Heart Study") +
  theme_minimal()
p5

# Save the plot
ggsave(p5, file = here("figures/", "Plot_A6.png"))

# Performing Fisher exact test
crosstab %>%
  table() %>%
  fisher.test()

# Cleaning up the environment
rm(crosstab)


## 3. How much do the following risk factors contribute to explaining occurrence
## of coronary heart disease? Are these contributions different between men and women?
## - total cholesterol, smoking status, age category, systolic blood pressure, BMI,
## and prevalence of hypertension
# Looping through the risk factors with Chi-square test
# total cholesterol, smoking status, agecategory, systolic blood pressure,
# BMI, and prevalence of hypertension
risk_factors <- list("totChol", "currentSmoker", "age_cat", "sysBP", "BMI",
                     "prevalentHyp")

# Creating empty list to store results
significance_tb <- list()

# Looping through the risk factors
for (risk_factor in risk_factors) {
  # Looping through both sexes and both gender
  for (sex_var in list("female", "male", c("female", "male"))) {
    # Creating crosstab for each risk factor
    crosstab <- framingham_df %>%
      mutate(age_cat = cut(age,
        breaks = c(39, 59, 79),
        labels = c("40-59", "60-79")
      )) %>%
      drop_na(all_of(risk_factor)) %>%
      select(sex, all_of(risk_factor), anyCHD)

    # Setting Chi-square test for categorical variables
    if (class(crosstab[[risk_factor]]) == "factor") {
      test_type <- "chisq.test"
    }

    # Tests for continuous variables
    else { # nolint
      # Testing for normal distribution
      norm_dist <- crosstab %>%
        pull(all_of(risk_factor)) %>%
        shapiro.test()

      # T-test for continuous variables normally distributed
      if (norm_dist$p.value > 0.01) {
        test_type <- "t.test"
      }

      # Wilcoxon test for continuous variables, not normally distributed
      else { # nolint
        test_type <- "wilcox.test"
      }
    }

    # Performing the test
    result <- suppressWarnings({
      crosstab %>%
        filter(sex %in% sex_var) %>%
        select(anyCHD, all_of(risk_factor)) %>%
        table() %>%
        get(test_type)()
    })

    # Cleaning up the environment
    rm(crosstab)

    # Outputting the result
    cat(paste0("Risk factor: ", risk_factor, ", Gender: ",
               paste(sex_var, collapse = " & ")))

    # Setting significance level to 0.01
    if (result$p.value < 0.01) {
      significance_tb <- append(significance_tb, list(c(risk_factor,
                                  paste(sex_var, collapse = " & "), "Yes"
                                )))
      cat(" -> Significant risk to attain a heart disease (", test_type, ")\n")
    } else {
      significance_tb <- append(significance_tb, list(c(risk_factor,
                                  paste(sex_var, collapse = " & "), "No"
                                )))
      cat(" -> No significant risk to attain a heart disease()", test_type, ")\n")
    }
  }
}

# Cleaning up the environment
rm(framingham_df, risk_factors, risk_factor, sex_var, test_type, result)

# Viewing the list
View(significance_tb)

# Transforming the list into a dataframe
significance_tb_df <- significance_tb %>%
  as.data.frame()

# Giving the columns numbers (1, 2, 3 etc.)
colnames(significance_tb_df) <- c(1 : length(significance_tb_df)) # nolint

# Giving the rows names
rownames(significance_tb_df) <- c("RiskFactor", "Sex", "Significant")

# Viewing the dataframe
View(significance_tb_df)

# Transforming the dataframe into tibble
significance_tb_tidied <- significance_tb_df %>%
  t() %>%
  as_tibble()

# Viewing the tibble
View(significance_tb_tidied)

# Creating table A1
t1 <- significance_tb_tidied %>%
  group_by(Significant) %>%
  mutate(Group = Significant) %>%
  gt() %>%
  tab_header(title = "Table A1: Risk factors",
             subtitle = "in relation to occurrence of coronary heart disease") %>%
  tab_source_note(source_note = md("*Data source: Framingham*")) %>%
  cols_hide(columns = Group) %>%
  cols_label(Group = md("*Significance*"), RiskFactor = md("**Risk Factor**"),
             Sex = md("**Gender**")) %>%
  cols_move(columns = c(RiskFactor, Sex), after = Group) %>%
  cols_align(align = "left", columns = everything()) %>%
  tab_spanner(label = md("*grouped by Significance*"), columns = RiskFactor) %>%
  tab_style(
    style = list(cell_fill(color = "#bd5c5c")),
    locations = cells_body(rows = Significant == "Yes")
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#4b8b4b")),
    locations = cells_body(rows = Significant == "No")
  )
t1
gtsave(t1, filename = "Table_A1.pdf", path = here("tables/"))


# Creating table A2
t2 <- significance_tb_tidied %>%
  group_by(Sex) %>%
  mutate(Group = Sex) %>%
  gt() %>%
  tab_header(title = "Table A2: Risk factors",
             subtitle = "in relation to occurrence of coronary heart disease") %>%
  tab_source_note(source_note = md("*Data source: Framingham*")) %>%
  cols_hide(columns = Group) %>%
  cols_label(Significant = md("**Significance**"), RiskFactor = md("**Risk Factor**")) %>%
  cols_move(columns = c(RiskFactor, Sex), after = Group) %>%
  cols_align(align = "left", columns = everything()) %>%
  tab_spanner(label = md("*grouped by Gender*"), columns = Significant) %>%
  tab_style(
    style = cell_fill(color = "#bd5c5c"),
    locations = cells_body(rows = Significant == "Yes")
  ) %>%
  tab_style(
    style = cell_fill(color = "#4b8b4b"),
    locations = cells_body(rows = Significant == "No")
  )
t2
gtsave(t2, filename = "Table_A2.pdf", path = here("tables/"))

# Cleaning up the environment
rm(significance_tb, significance_tbsignificance_tb_df, significance_tb_tidied, t1, t2)
