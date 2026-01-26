# Modeling River-floodplain exchange

This repo contains code to run our machine learning analysis of United States river-floodplain exchange.

## Setup

Use `environment.yaml` to generate conda environment.

You will need a dev version (at the time) of `terra` and `flowdem`. We used `terra` version 1.8-81 and `flowdem` version 0.2. Download within an interactive R session within the conda environment via:

```r
remotes::install_github("rspatial/terra",  repos='https://rspatial.r-universe.dev')
remotes::install_github("KennethTM/flowdem")
```

The data repo is specified by a simulink at `~/data/path_to_data`. Other lookup tables auxillary datasets are stored directly in `~/data`. Make sure you inspect the functions to understand how the NHD-HR hydrography, DEMs, and other datasets must be stored to run. Note that setting up your own data repo is non-trivial.

## To run

Make sure you are familiar with the `targets` R pipelining tool and `crew` distributed computing helpers. Realistically, this workflow requires the use of an HPC server and will launch over 60,000 targets.

`run.sh` will launch a non-interactive R session and run the pipeline specifications in `~/src/runTargets.R`. It is written specifically for our SLURM scheduler- your HPC will require different resource settings.

This repo is setup to run in parallel across processors on a single machine. Parallel computing settings are specified in the targets pipeline (`_targets.R`). You will need to manually adjust the parallel computing settings for each task's resource needs and given the constraints of your HPC. This requires familiarity with the pipeline, the resource requirements of each task, and your HPC's available processors. A fully distributed workflow across machines is possible using `crew`, but we avoided it to limit the number of jobs we submit at a time.