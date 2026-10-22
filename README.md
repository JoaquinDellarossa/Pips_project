# PowerLab

PowerLab is an interactive Shiny app for learning and planning statistical power analyses.

Before collecting data, psychologists need to decide how many participants to include; with too few, a real effect can easily be missed. PowerLab lets students and researchers explore how sample size, effect size, and significance level determine statistical power, and calculates the sample size needed to reach a target power. Power is calculated analytically from the noncentral *t* and *F* distributions for:

- One-sample *t*-test
- Independent two-sample *t*-test
- One-way ANOVA
- Two-way ANOVA (main effect A, main effect B, or the A × B interaction)

## Getting started

Clone the repository:

```bash
git clone https://github.com/JoaquinDellarossa/Pips_project.git
cd Pips_project
```

PowerLab requires R 4.4 or newer; RStudio is recommended. Install the packages needed to open an R Markdown document:

```r
install.packages(c("knitr", "rmarkdown"))
```

Open `main.rmd` in RStudio and run its chunks from top to bottom. The first chunk installs any missing app packages. The app chunk is not evaluated while knitting and should be run manually.

## Repository structure

```
Pips_project/
├── src/
│   └── power_lab_app.R  # Shiny app: power calculations, user interface, and plots
├── main.rmd             # Main deliverable: installs packages and runs the app
└── README.md
```

## Using the project

`main.rmd` explains the project and starts the app. The app can also be started outside the document:

```bash
Rscript -e "shiny::runApp('src/power_lab_app.R')"
```

Inside the app:

1. Follow the short introduction (why power matters, key concepts, the four tests, and an interactive preview), or skip it.
2. Choose a learning level: *Beginner* shows the results and power curve; *Advanced* also shows degrees of freedom, the noncentrality parameter, critical values, and the sampling distributions. You can switch at any time.
3. On the Explore page, choose a test, set the sample size, effect size (Cohen's *d* or *f*), α, and target power, and read off the power, total sample size, and required *n*.

Every ⓘ button and most chart elements open a short explanation.

## Dependencies

The app uses `shiny` for the interface and `ggplot2` (version 3.4.0 or newer) for the plots. `knitr` and `rmarkdown` support the main document. `main.rmd` installs `shiny` and `ggplot2` when they are missing.

## Reference

Cohen, J. (1988). *Statistical power analysis for the behavioral sciences* (2nd ed.). Lawrence Erlbaum Associates.
