###################################################################################################
###################################################################################################
## POXVIRUS MANUSCRIPT ANALYSES
###################################################################################################
###################################################################################################


###################################################################################################
## PACKAGES
###################################################################################################

library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(scales)
library(ggpubr)
library(arm)


###################################################################################################
###################################################################################################
## 1. PREVALENCE MODEL
###################################################################################################
###################################################################################################

## Load final dataset used for prevalence analysis
## Replace with the final file you will provide
# final_prevs_tot_subset <- read.csv("/Users/IOC/Downloads/Final.csv")


## Make sure variables are correctly formatted
final_prevs_tot_subset$BeforeAfterQuarantine <- factor(
  final_prevs_tot_subset$BeforeAfterQuarantine,
  levels = c("Before", "After")
)

final_prevs_tot_subset$Treatment <- factor(
  final_prevs_tot_subset$Treatment,
  levels = c("Urban", "Suburban", "Rural")
)



###################################################################################################
## Weighted linear model
###################################################################################################

modelo_prev <- lm(
  sqrt(Prevalence) ~ Treatment + BeforeAfterQuarantine,
  weights = Quantidade,
  data = final_prevs_tot_subset
)

modelo_prev_null <- lm(
  sqrt(Prevalence) ~ 1,
  weights = Quantidade,
  data = final_prevs_tot_subset
)


## Model comparison
anova(modelo_prev, modelo_prev_null)


## Model results
summary(modelo_prev)
drop1(modelo_prev, test = "F")
confint(modelo_prev)


###################################################################################################
## Prevalence boxplot: Before vs After
###################################################################################################

p_prevalence <- ggplot(
  final_prevs_tot_subset,
  aes(
    x = BeforeAfterQuarantine,
    y = Prevalence,
    fill = BeforeAfterQuarantine
  )
) +
  geom_boxplot(
    width = 0.55,
    alpha = 0.75,
    outlier.size = 2
  ) +
  scale_fill_manual(
    values = c(
      "Before" = "#F2B66D",
      "After" = "#6BA6CC"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    x = "COVID-19 Anthropause",
    y = "Poxvirus prevalence"
  ) +
  theme_classic(base_size = 15) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black", size = 13),
    legend.position = "none"
  )

p_prevalence


###################################################################################################
###################################################################################################
## 2. INDIVIDUAL POXVIRUS OCCURRENCE MODEL
###################################################################################################
###################################################################################################

## Assumes data_final_subset is the final individual-level dataset

data_final_subset$BeforeAfterQuarentine <- factor(
  data_final_subset$BeforeAfterQuarentine,
  levels = c("Before", "After")
)

data_final_subset$Site2 <- factor(
  data_final_subset$Site2,
  levels = c("Urban", "Suburban", "Rural")
)


###################################################################################################
## Descriptive prevalence by habitat
###################################################################################################

prevalence_by_site <- data_final_subset %>%
  group_by(Site2) %>%
  summarise(
    n_total = n(),
    n_infected = sum(PoxPresence == 1, na.rm = TRUE),
    prevalence = n_infected / n_total * 100,
    .groups = "drop"
  )

print(prevalence_by_site)


###################################################################################################
## Binomial model
###################################################################################################

modelo <- glm(
  PoxPresence ~ Site2 + BeforeAfterQuarentine,
  data = data_final_subset,
  family = binomial
)

modelo_null <- glm(
  PoxPresence ~ 1,
  data = data_final_subset,
  family = binomial
)


## Model comparison
anova(modelo, modelo_null, test = "Chisq")
AIC(modelo, modelo_null)


## Model results
summary(modelo)
drop1(modelo, test = "Chisq")
confint(modelo)

###################################################################################################
## Individual occurrence plot: Before vs After
###################################################################################################

p_occurrence <- ggplot(
  data_final_subset,
  aes(
    x = BeforeAfterQuarentine,
    fill = PoxPresence
  )
) +
  geom_bar(
    position = "fill",
    width = 0.62
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.03))
  ) +
  scale_fill_manual(
    values = c(
      "Uninfected" = "#4C78A8",
      "Infected" = "#E45756"
    )
  ) +
  labs(
    x = "COVID-19 Anthropause",
    y = "Proportion of individuals",
    fill = "Poxvirus occurrence"
  ) +
  theme_classic(base_size = 15) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black", size = 13),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12)
  )

p_occurrence


###################################################################################################
## Combined Figure 2
## A = occurrence
## B = prevalence
###################################################################################################

figure2 <- ggarrange(
  p_occurrence,
  p_prevalence,
  ncol = 1,
  nrow = 2,
  labels = c("A", "B"),
  font.label = list(
    size = 18,
    face = "bold"
  )
)

figure2


###################################################################################################
###################################################################################################
## 3. 2020 SAMPLING SENSITIVITY
###################################################################################################
###################################################################################################

pox_year <- data_final_subset %>%
  filter(
    !is.na(Year),
    !is.na(PoxPresence)
  ) %>%
  group_by(Year) %>%
  summarise(
    n = n(),
    positive = sum(PoxPresence == 1),
    prevalence = positive / n,
    .groups = "drop"
  )

print(pox_year)


## 2020 sample
pox_2020 <- pox_year %>%
  filter(Year == 2020)

n_2020 <- pox_2020$n
positive_2020 <- pox_2020$positive


## Exact 95% binomial CI
binom_2020 <- binom.test(
  x = positive_2020,
  n = n_2020
)

print(binom_2020)


## 2019 prevalence
prev_2019 <- pox_year %>%
  filter(Year == 2019) %>%
  pull(prevalence)


## Pooled 2017–2019 prevalence
pre_covid <- data_final_subset %>%
  filter(
    Year %in% 2017:2019,
    !is.na(PoxPresence)
  ) %>%
  summarise(
    n = n(),
    positive = sum(PoxPresence == 1),
    prevalence = positive / n
  )

print(pre_covid)


## Detection probabilities
prob_zero_2019 <- dbinom(
  0,
  size = n_2020,
  prob = prev_2019
)

prob_detect_2019 <- 1 - prob_zero_2019


prob_zero_pre_covid <- dbinom(
  0,
  size = n_2020,
  prob = pre_covid$prevalence
)

prob_detect_pre_covid <- 1 - prob_zero_pre_covid


## Summary
sensitivity_results <- tibble(
  Reference = c(
    "2019 prevalence",
    "2017–2019 pooled prevalence"
  ),
  Assumed_prevalence = c(
    prev_2019,
    pre_covid$prevalence
  ),
  N_2020 = n_2020,
  Prob_zero_positives = c(
    prob_zero_2019,
    prob_zero_pre_covid
  ),
  Prob_detect_at_least_one = c(
    prob_detect_2019,
    prob_detect_pre_covid
  )
)

print(sensitivity_results)


###################################################################################################
###################################################################################################
## 4. KAPLAN–MEIER ANALYSIS AFTER THE ANTHROPAUSE
###################################################################################################

pox_after <- data_final2 %>%
  filter(BeforeAfterQuarentine == "After")

pox_after$DaysSinceAnthropause <- as.numeric(
  as.Date(pox_after$Date2) - as.Date("2020-06-01")
)

## Survival analysis
# Create survival object
pox_1yr <- subset(pox_after, DaysSinceAnthropause <= 500)
pox_1yr$event <- pox_1yr$PoxPresence == 1 
surv_obj <- Surv(time = pox_1yr$DaysSinceAnthropause, event = pox_1yr$event)

# Fit survival model
surv_obj <- Surv(time = pox_1yr$DaysSinceAnthropause, event = pox_1yr$event)
km_fit   <- survfit(surv_obj ~ Site2, data = pox_1yr)
summary(km_fit)
# Prepare survival summary
ss <- survminer::surv_summary(km_fit, data = pox_1yr) %>%
  mutate(
    cuminc = 1 - surv,
    cil    = 1 - upper,
    ciu    = 1 - lower
  )

# Clean legend labels (remove 'Site2=')
ss$strata <- gsub("Site2=", "", ss$strata)

# Find maximum cumulative incidence
ymax <- max(ss$cuminc, na.rm = TRUE)

ggplot(ss, aes(x = time, y = cuminc, color = strata, fill = strata)) +
  geom_ribbon(aes(ymin = cil, ymax = ciu), alpha = 0.25, color = NA) +
  geom_step(linewidth = 1) +
  coord_cartesian(
    xlim = c(200, max(ss$time, na.rm = TRUE)),
    ylim = c(0, ymax)
  ) +
  labs(
    x = "Days Since Anthropause",
    y = "Cumulative pox incidence (1 - S(t))",
    color = "Habitat",
    fill = "Habitat"
  ) +
  scale_color_manual(values = c("red", "green3", "blue")) +
  scale_fill_manual(values = c("red", "green3", "blue")) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 20, 10, 10),
    legend.title = element_blank(),                # removes legend title
    legend.position = c(0.8, 0.7),                 # move legend inside plot
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_blank()
  )

## End of code