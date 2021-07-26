## Introduction

The repository contains the scripts used in preparing the **Rates of contributory de novo mutation in high and low-risk autism families** by Seungtai Yoon, Adriana Munoz, et al. manuscript, soon to be published in _Nature Communications Biology_.

## Description of the contents

The **input** directory contains copies of the Supplementary Data associated with the manuscript. A detailed description of these data files can be obtained from _Nature Communications Biology_'s copy of the manuscript. The **src** directory python scripts use these data files and perform all the analyses presented in the manuscript. As a result, the scripts generate the tables and figures shown in the manuscript and the associated supporting information. As a convenience, the repository also stores the **results** of the execution of the scripts. Thus, readers can examine the scripts together with their inputs and outputs without setting up an execution environment for themselves. 

We also provide instructions for executing the scripts for the readers who want to run the scripts to re-produce our analysis or apply similar analysis over their data. 

## Instructions for running the analysis

#### Setting up the environment

We wrote our python scripts in Python 2.7 using few packages outside of the standard python library. You must install precisely the versions of the external packages indicated in the **environment.yml** file to ensure that the scripts generate output identical to the one stored in the **results** directory. There are numerous ways you can achieve that. Here we will provide detailed instructions using the conda package manager, the one we use to manage our execution environments.

First, you have to install conda following the instructions in [Conda Installation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

Then you should create and activate a conda environment using the specification in the **environment.yml" file. You can do that with the following commands executed from the directory containing the clone of this repository:

    $ cd .../denovoInHighAndLowRiskPaper
    $ conda env create
    $ conda activate denovoInHighAndLowRiskPaper

**NOTE**: The **denovoInHighAndLowRiskPaper** name of the conda environment is specified in the **environment.yml** file. 

#### Execute scripts one by one

After you activate the **denovoInHighAndLowRiskPaper** environment, you would be able to run each of the scripts using a command like:

    $ python src/figRatesVsAge.py

#### Execute all the scripts using **make**

The **Makefile** contains the instructions on how to execute the scripts needed to generate all the results. The standard tool **make** can read and execute the instructions included in the **Makefile**. **make** is installed by default on Linux distributions, is very easy to install on Mac, and, if needed, is available on Windows (i.e., http://gnuwin32.sourceforge.net/packages./make.htm)

The following simple commands would regenerate all the results in a new directory called **myResults**. 

    $ mkdir myResults
    $ cd myResults
    $ make -f ../Makefile -j

You can then compare the newly generated files in **myResults** to those stored in the existing **results** directory: if the inputs and the scripts have not changed, the two sets of result files should be identical. 

## A guide to the scripts

#### Table scripts

We generate **Tables 1-5** and **Supplementary Table 1** in two steps. The first step is the heavy one, and we implemented it in the **result_tables.py** script. It computes all the necessary results and stores them in four intermediate files. We implemented the second step using five additional python scripts that generate the six tables included in the manuscript with their precise content and layout. 

| Manuscript table      | Intermediate file                   | Second step script                             |
|-----------------------|-------------------------------------|------------------------------------------------|
| Table 1               | small_scale_result_table.txt        | tabSmallScale.py                               |
| Table 2               | CNV_result_table.txt                | tabAllCNVs.py                                  |
| Table 3               | mered_simplex_multiplex_results.txt | tabSmallScaleAndCNVs.py                        |
| Table 4               | CNV_result_table.txt                | tabOneGeneCNVs.py                              |
| Table 5               | intronic_result_table.txt           | tabIntronicPeripheral.py inter-coding_intronic |
| Supplementary Table 1 | intronic_result_table.txt           | tabIntronicPeripheral.py peripheral            |

sReal.Cu
sReal.Su
sReal.Nu
sReal.Ca
sReal.Sa
sReal.Na
sReal.RSa
sReal.ESa
sReal.delta
sReal.AD
bcAD.left95
bcAD.right95
esAD.pvOne
esAD.z
esAD.pvOneAn
sReal.PC
bcPC.left95
bcPC.right95


#### Figure scripts

**Figure 1** and **Supplementary Figures 1-7** are each generated by a dedicated script:

| Figure                 | script                              |
|------------------------|-------------------------------------|
| Figure 1               | figChildrenScatter.py               |
| Supplementary Figure 1 | figCNVsByFilter.py                  |
| Supplementary Figure 2 | figSnvPower.py                      |
| Supplementary Figure 3 | figParentalAges.py                  |
| Supplementary Figure 4 | figCNVPower.py                      |
| Supplementary Figure 5 | figPCbyCNVgenes.py <sup>1</sup>     |
| Supplementary Figure 6 | figRatesVsAge.py                    |
| Supplementary Figure 7 | figISBSignalPower.py <sup>1,2</sup> |

You can execute each of the figure scripts with or without an argument. If you provide one, the script will use it as the file's name to store the figure. The extension of the file will determine the image format. Allowed extensions include ".pdf", ".png", and ".jpeg". If you provide no argument, the script will create a 'live' figure that you can interact with (i.e., you can resize, zoom in or out, save the figure, etc.).

**NOTES:** 

<sup>1</sup> The script uses a random number generator and will generate slightly different results in every execution. If stable results are required, you may pass a number used to initialize the random number generator as a second argument to the script. 

<sup>2</sup> The script uses information from the **intronic_result_table.txt** generated as part of the table generation process. **figISBSignalPower.py** will look for the **intronic_result_table.txt** in the current working directory. You can generate this file by running the **result_tables.py** script. Alternatively, you can switch your current working directory to the **results** directory that already contains the **intronic_result_table.txt** file.

#### Variant property analysis

Finally, **Supplementary Figure 8-32** and **Supplementary Table 2** are all generated by the **draw_float_property_hist.py** script. 