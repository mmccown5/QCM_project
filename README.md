# QCM_project

For QCM Class, we will complete task1 of the NeuroIPS competition. Full documentation for the competition can be found online 
at [openproblems.bio/neurips_docs/](https://openproblems.bio/neurips_docs/). This repository began with the starter kit found on the [Quickstart](https://openproblems.bio/neurips_docs/submission/quickstart/) page.

## Method
To be determined

## Folder Structure
​
```
├── LICENSE                                 # MIT License
├── README.md                               # This repo outline
├── starter_kit_README.md                   # Starter information from the starter kit
├── .gitignore                              # 
├── bin/                                    # Binaries needed to generate a submission
│   ├── check_format
│   ├── nextflow
│   └── viash
├── config.vsh.yaml                         # Viash configuration file
├── script.py                               # Script containing your method
├── sample_data/                            # Small sample datasets for unit testing and debugging
│   ├── openproblems_bmmc_cite_starter/     # Contains H5AD files for CITE data
│   └── openproblems_bmmc_multiome_starter/ # Contains H5AD files for multiome data
├── benchmark_data/                         # Benchmark or exploratory datasets
│   ├── cite/                               # Contains H5AD files for CITE data
│   ├── multiome/                           # Contains H5AD files for multiome data
│   ├── LICENSE                             # MIT License
│   └── README.md                           # See https://openproblems.bio/neurips_docs/data
├── scripts/                                # Scripts to test, generate, and evaluate a submission
│   ├── 0_sys_checks.sh                     # Checks that necessary software installed
│   ├── 1_unit_test.sh                      # Runs the unit tests in test.py
│   ├── 2_generate_submission.sh            # Generates a submission pkg by running your method on validation data
│   ├── 3_evaluate_submission.sh            # (Optional) Scores your method locally
│   └── nextflow.config                     # Configurations for running Nextflow locally
└── test.py                                 # Default unit tests. Feel free to add more tests, but don't remove any.
```
