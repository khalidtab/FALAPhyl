conda activate $CONDA_PREFIX
R -e 'if (!requireNamespace("devtools")) install.packages("devtools") ; devtools::install_github("Russel88/DAtest")'