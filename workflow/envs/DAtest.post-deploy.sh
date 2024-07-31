conda activate $CONDA_PREFIX
R -e 'if (!requireNamespace("devtools")) install.packages("devtools") ; devtools::install_github("Russel88/DAtest", upgrade = "never")'
R -e 'devtools::install_version("statmod", version = "1.4.33", repos ="http://cran.us.r-project.org”, upgrade = "never") ; devtools::install_version("cplm", version = "0.7-5", repos = "http://cran.us.r-project.org", upgrade = "never")'
R -e 'devtools::install_github("himelmallick/Tweedieverse", upgrade = "never")'