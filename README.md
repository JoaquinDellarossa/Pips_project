# Pips_project

This repository contains PowerLab, an interactive Shiny app for learning about statistical power and planning sample sizes. We created it for our group project in the Programming in Psychological Science (PIPS) course of the Research Master Psychology at the University of Amsterdam.

PowerLab lets users calculate prospective power for a one-sample t-test, an independent two-sample t-test, a one-way ANOVA, and a two-way between-subjects ANOVA. Users can see how sample size, effect size (Cohen's *d* or *f*) and significance level shape the sampling distributions under H₀ and H₁, follow power along a power curve, and find the sample size needed to reach a target power.

## Citation

[Add the citation for your repository here.]

## Getting started

Clone the repository:

```bash
git clone https://github.com/JoaquinDellarossa/Pips_project.git
cd Pips_project
```

The app requires R; RStudio is recommended. Install the required packages:

```r
install.packages(c("shiny", "ggplot2"))
```

Open RStudio in the `Pips_project` folder and start the app from the console:

```r
shiny::runApp("Power_Lab_with_ANOVA.R")
```

The app opens in your browser or in the RStudio viewer. To stop it, press Esc in the console or close the window.

## Repository structure

```
Pips_project/
├── Power_Lab_with_ANOVA.R   # PowerLab Shiny app
├── Power Lab.R              # [Describe this file]
├── Pips_project.R           # [Describe this file]
└── README.md
```

## Using the app

When PowerLab opens, a short four-step introduction explains why power matters and lets you choose a learning level. Choose **Skip introduction** to go straight to the lab.

- **Beginner** shows the core inputs, the power result, the sampling distributions and the power curve.
- **Intermediate** adds the degrees of freedom, noncentrality parameter and critical values, a plot of the α, β and power regions, and a sample-size comparison.
- **Expert** adds a sensitivity analysis and downloads of an analysis report and a power-curve CSV.

You can change the level at any time; the calculations stay the same.

In the lab, choose a test, set the sample size, effect size, α and target power, and read off the power, total N and required n. For the ANOVAs, *n* is the sample size per group (one-way) or per cell (two-way). Click any ⓘ button, or click directly on a curve, boundary or shaded area in the plots, for an explanation of what it represents.

## Dependencies

PowerLab uses Shiny for the interactive interface and ggplot2 (version 3.4.0 or newer) for the plots. Power is computed from the noncentral *t* and *F* distributions in base R.

## Reference

Cohen, J. (1988). *Statistical power analysis for the behavioral sciences* (2nd ed.). Lawrence Erlbaum Associates.
