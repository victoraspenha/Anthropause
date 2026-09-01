# Anthropause
Disappearance of a common wildlife disease during the COVID-19 pandemic, followed by its return first in human-developed areas   

POXVIRUS ANTHROPAUSE ANALYSIS
=============================

This R script contains the statistical analyses and figures used in the manuscript examining changes in avian poxvirus infection associated with the COVID-19 Anthropause.

The script has been cleaned to retain only analyses relevant to the manuscript. Exploratory analyses, maps, forest plots, alternative models, and unused figures have been removed.


REQUIRED R PACKAGES
-------------------

dplyr
ggplot2
survival
survminer
scales
ggpubr
arm


DATASETS
--------

1. final_prevs_tot_subset

Sampling-event-level dataset used for the poxvirus prevalence analysis.

Main variables:

- Date: sampling date
- Prevalence: proportion of infected individuals during a sampling event
- Quantidade: number of individuals sampled
- Treatment: habitat category (Urban, Suburban, Rural)
- BeforeAfterQuarantine: sampling period relative to the COVID-19 Anthropause (Before, After)


2. data_final_subset

Individual-level dataset used for occurrence, sensitivity, and survival analyses.

Main variables:

- BandID: individual bird identifier
- Date2: sampling date
- Year: sampling year
- Site2: habitat category (Urban, Suburban, Rural)
- PoxPresence: poxvirus infection status (Infected, Uninfected)
- BeforeAfterQuarentine: sampling period relative to the COVID-19 Anthropause (Before, After)

Note: "BeforeAfterQuarentine" retains the original spelling used in the individual-level dataset, whereas the prevalence dataset uses "BeforeAfterQuarantine".


ANALYSES
--------

1. POXVIRUS PREVALENCE

Response:
Square-root-transformed poxvirus prevalence.

Model:
Weighted linear model.

Formula:

sqrt(Prevalence) ~ Treatment + BeforeAfterQuarantine

Sampling events are weighted by the number of birds sampled (Quantidade).

The model evaluates whether poxvirus prevalence differs:

- among Urban, Suburban, and Rural habitats
- before versus after the COVID-19 Anthropause

The script includes:

- histogram of square-root-transformed prevalence
- comparison with an intercept-only model
- model summary
- drop1 tests
- 95% confidence intervals
- model diagnostic plot


2. POXVIRUS OCCURRENCE

Response:
Individual poxvirus infection status.

Model:
Binomial generalized linear model.

Formula:

PoxPresence ~ Site2 + BeforeAfterQuarentine

The model evaluates whether individual infection probability differs:

- among Urban, Suburban, and Rural habitats
- before versus after the COVID-19 Anthropause

The script includes:

- descriptive prevalence by habitat
- comparison with an intercept-only model
- model summary
- likelihood-ratio tests
- 95% confidence intervals
- model diagnostic plot


3. 2020 SAMPLING SENSITIVITY ANALYSIS

This analysis evaluates whether the absence of detected poxvirus infections in 2020 could be explained by reduced sampling effort.

The analysis calculates:

- annual sample sizes
- annual numbers of infected birds
- annual prevalence
- exact 95% binomial confidence interval for 2020 prevalence
- probability of detecting at least one infected bird in 2020 if prevalence had remained at the 2019 level
- probability of detecting at least one infected bird in 2020 if prevalence had remained at the pooled 2017–2019 level

Detection probability is calculated as:

1 - (1 - p)^n

where:

p = reference prevalence
n = number of birds sampled in 2020


4. KAPLAN-MEIER / CUMULATIVE INCIDENCE ANALYSIS

This analysis examines the timing of poxvirus infections following the Anthropause.

Only observations during the first 500 days after June 1, 2020 are included.

Time variable:

DaysSinceAnthropause =
sampling date - June 1, 2020

Event definition:

1 = Infected
0 = Uninfected

Kaplan-Meier curves are estimated separately for:

- Urban
- Suburban
- Rural

Survival estimates are converted to cumulative incidence:

Cumulative incidence = 1 - S(t)

A log-rank test is used to compare the cumulative incidence curves among habitats.

The manuscript figure shows cumulative poxvirus incidence over time with 95% confidence intervals.


FIGURES
-------

FIGURE 2

Panel A:
Proportion of individual birds classified as infected or uninfected before and after the Anthropause.

Panel B:
Distribution of sampling-event-level poxvirus prevalence before and after the Anthropause.

Panels are arranged vertically to improve readability.


CUMULATIVE INCIDENCE FIGURE

Kaplan-Meier estimates transformed to cumulative incidence (1 - S(t)) for Urban, Suburban, and Rural habitats.

Shaded areas represent 95% confidence intervals.


IMPORTANT CODING NOTES
----------------------

PoxPresence is stored as:

"Infected"
"Uninfected"

For the survival analysis it is converted to:

Infected   = 1
Uninfected = 0

Before/after variable names differ between datasets:

Individual-level data:
BeforeAfterQuarentine

Prevalence data:
BeforeAfterQuarantine

These names should not be interchanged.
