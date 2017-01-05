PHRAPL --  user manual
=======

##**WARNING**
This user manual is under development. Please send me an [email](ariadna.biologia@gmail.com) if you have comments.
Figures and examples need to be added.

________

<img src="https://github.com/ariadnamorales/phrapl-manual/blob/master/phrapl_logo.png?raw=true" width="200" height="200" />

### PHRAPL project [web site](http://www.phrapl.org/)

### CITATION
- Jackson N, Morales AE, Carstens BC, O'Meara BC (2016) PHRAPL: Phylogeographic Inference using Approximate likelihoods. Systematic Biology. In press

- Jackson N, Carstens BC, Morales AE, O’Meara BC (2016) Species delimitation with gene flow. Systematic Biology. [DOI: 10.1093/sysbio/syw117](http://sysbio.oxfordjournals.org/content/early/2016/12/05/sysbio.syw117.short?rss=1)

- Morales AE, Jackson N, Dewey T, O’Meara BC, Carstens BC (2016) Speciation with gene flow in North American *Myotis* bats. Systematic Biology. [DOI: 10.1093/sysbio/syw100](http://sysbio.oxfordjournals.org/content/early/2016/12/20/sysbio.syw100.full?sid=4376fa06-5da3-4e21-9990-faf3c57b9298)

### CODE
`PHRAPL` is written in `R`, but it uses `perl` and `ms` to perform simulations. The pre-CRAN (code under development) can be found in [github.](https://github.com/bomeara/phrapl)

## Why to use `PHRAPL`?
Phylogeographic research aims to understand the recent history of species. Over the last decades, researchers have increasingly incorporated demographic models in order to estimate parameters (i.e., divergence times, population sizes, and rates of migration and expansion) that can contribute phylogeographic inference. 
Typically, this is conducted via the use of software packages that contain specified models (n-island models or fixed topologies). Alternatively, simulation-based approaches allow researchers to customize models for the particular details of their system, and may be useful in testing preexisting biogeographic hypotheses. Because the demographic model is central to the analysis in either case, researchers may wish to assess the appropriateness of their model to the data. `PHRAPL` is designed to give such a tool to researchers. 

## How `PHRAPL` works?
`PHRAPL` simulates genealogies under a wide range of demographic models and compares the empirical genealogies to the simulated gene tree distributions. Demographic models that are probable given the data will contain many genealogies that match the estimated gene trees. Because the proportion of matching gene trees for a given model is equivalent to the probability of the data given the model and parameter values, we can use this value in an information theoretic framework to evaluate the relative weight of all models. This provides the researcher with an independent assessment of both the best model, given the data as well as the ability to calculate the model likelihoods of classes of models (e.g., n-island vs. isolation models).

Watch [these YouTube videos](https://www.youtube.com/watch?v=UC4Mj1K6c0k) to learn more about `PHRAPL`.

## Contents

-    [CRAN-version Vignette](https://github.com/bomeara/phrapl/blob/master/doc/phrapl_vignette.Rmd)

1. [Installation](https://github.com/ariadnamorales/phrapl-manual/blob/master/1.Installation.Rmd)
  - [Phydocker for windows users](http://www.phrapl.org/#distribution): It lets you use `PHRAPL` without having to install R, perl, ms, etc, just docker. Recommended only for testing.
2. [Input files](https://github.com/ariadnamorales/phrapl-manual/blob/master/2.Input_files.Rmd)
3. [Generate a set of models (migrationArray object)](https://github.com/ariadnamorales/phrapl-manual/blob/master/3.Generate_set_of_models.Rmd)
  - [3a. How models are built?](https://github.com/ariadnamorales/phrapl-manual/blob/master/3a.How_models_are_built.Rmd)
  - [3b. Plotting models](https://github.com/ariadnamorales/phrapl-manual/blob/master/3b.Plotting_models.Rmd)
4. [Subsampling and creating an input for `PHRAPL`](https://github.com/ariadnamorales/phrapl-manual/blob/master/4.Subsample_CreateInput.Rmd)
5. [Running `PHRAPL` (GridSearch option)](https://github.com/ariadnamorales/phrapl-manual/blob/master/5.Run_Phrapl.Rmd)
6. [Post-processing](https://github.com/ariadnamorales/phrapl-manual/blob/master/6.Post-processing.Rmd)

Advanced topics:
7. [Sensitivity Analyses](https://github.com/ariadnamorales/phrapl-manual/blob/master/7.SensitivityAnalyses.Rmd)
8. [Species delimitation]()


### Do you have a question about `PHRAPL` or want to report a bug?
Post it in the [phrapl-users](https://groups.google.com/forum/#!forum/phrapl-users) google group.