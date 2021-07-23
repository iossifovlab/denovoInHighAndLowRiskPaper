## Introduction

The repository contains the scripts used in preparing the **Rates of contributory de novo mutation in high and low-risk autism families** by Seungtai Yoon, Adriana Munoz, et al. manuscript, soon to be published in _Nature Communications Biology_.

## Description of the contents

The **input** directory contains copies of the Supplementary Data associated with the manuscript. A detailed description of these data files can be obtained from _Nature Communications Biology_'s copy of the manuscript. The **src** directory python scripts use these data files and perform all the analyses presented in the manuscript. As a result, the scripts generate the tables and figures shown in the manuscript and the associated supporting information. As a convenience, the repository also stores the **results** of the execution of the scripts. Thus, readers can examine the scripts together with their inputs and outputs without setting up an execution environment for themselves. 

We also provide instructions for executing the scripts for the readers who want to run the scripts to re-produce our analysis or apply similar analysis over their data. 

## Instructions for running the analysis

### Step 1. Setting up the environment

We wrote our python scripts in Python 2.7 using few packages outside of the standard python library. You must install precisely the versions of the external packages indicated in the **environment.yml** file to ensure that the scripts generate output identical to the one stored in the **results** directory. There are numerous ways you can achieve that. Here we will provide detailed instructions using the conda package manager, the one we use to manage our execution environments.

First, you have to install conda following the instructions in [Conda Installation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

Then you should create and activate a conda environment using the specification in the **environment.yml" files. You can do that with the following commands executed from the directory containing the clone of this repository:

    $ cd .../denovoInHighAndLowRiskPaper
    $ conda env create
    $ conda activate denovoInHighAndLowRiskPaper

### Step 2. Execute all the scripts using **make**
