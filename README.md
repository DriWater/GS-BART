# GS-BART: Bayesian Additive Regression Trees with Graph-split Decision Rules for Generalized Spatial Nonparametric Regressions

## Instructions for use

All workflow information to reproduce key results of the paper is contained in `reproduce.Rmd`. The steps for running the workflow are:

**Step 1**: Install all required R packages in `dependencies.txt`.

**Step 2** (Optional): Compile the source code for the GS-BART library by running the following command in your R Console:
```
install.packages('R packages/GSBart_0.0.1.tar.gz', type = "source", repo = NULL) 
install.packages('R packages/G2SBart_0.0.1.tar.gz', type = "source", repo = NULL)
```
This step is optional and only required if you want to re-run GS-BART model fitting and prediction.

**Step 3** (Optional): Modify the options in Line 19-25 of `reproduce.Rmd` accordingly. By default, key figures and tables are reproduced from pre-computed model outputs in `data/`. To re-generate model outputs, make the following change(s):
* To re-run GS-BART model in continuous simulation data, set `sim_continuous_fit = TRUE`.
* To re-run GS-BART model in count simulation data, set `sim_count_fit = TRUE`.
* To re-run GS-BART model in classification simulation data, set `sim_classification_fit = TRUE`.
* To re-run GS-BART model in continuous real data, set `real_continuous_fit = TRUE`.
* To re-run GS-BART model in count real data, set `real_count_fit = TRUE`.
* To re-run GS-BART model in classification real data, set `real_classification_fit = TRUE`.

**Step 4**: Run the workflow in `reproduce.Rmd` by compiling the R Markdown file. For example, this can be done by running the following command in your terminal:
```
Rscript -e "rmarkdown::render('reproduce.Rmd')"
```
Reproducible results can be found in the generated `reproduce.html` file.
