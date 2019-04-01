---
layout: page
title: Setup
permalink: /setup/
---

# Option 1: Launch Binder

For a zero installation option, launch the cell factory design course by clicking on the Binder button:
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/AgileBioFoundry/2019-04.17-cell-factory-design-brownbag/master?filepath=01-Intro-to-FBA.ipynb)

# Option 2: Install on your own machine

Below are the steps to install the cell factory design course on your own machine

## Prerequisites

Before installing any of the course-specific packages below, please head over to <https://pnnl-compbio.github.io/2019-03.14-15-PNNL-SWC/#python> and
follow the installation instructions for Python.

## Create the `cell-factory-design-course` conda environment

Please open a shell (terminal) and create a conda environment for the course using the following command.

    conda env create  -f environment.yml

If you're on Linux or OS X or are using the git-bash shell on Windows run

    source activate cell-factory-design-course

to activate the environment. If you're on Windows using the default terminal (cmd) run the following command instead.

    activate cell-factory-design-course

If you successfully activated you're environment, your command prompt will look similar to this.

    (cell-factory-design-course)$

Run the following command to check that the installation was successful (this can take a few seconds).

    python -c "from cameo import models;print(models.bigg.e_coli_core.optimize().objective_value)"

The output should be approximately `0.8739`.

Lastly, you need to add this virtual environment to the Jupyter notebook kernel

    python -m ipykernel install --user --name cell-factory-design-course --display-name "cell factory design course"

You also have received an email with with a download link for CPLEX (only intended for academic use). Download the respective installer for your platform and install CPLEX. After you succeeded follow the respective instructions for your platform.

### OS X and Linux

On Linux or OS X, run the following in your terminal (make sure you activated the `cell-factory-course` environment first) to install CPLEX for Python.

    pip install <path-to-your-cplex-installation>/CPLEX_Studio128/cplex/python/3.6/<platform>/

 For OS X, this should look like

 	pip install /Applications/CPLEX_Studio128/cplex/python/3.6/x86-64_osx/

### Windows

If you installed CPLEX on Windows, `cd` into the following directory first
    
    cd "C:\Program Files\IBM\ILOG\CPLEX_Studio128\cplex\python\3.6\<platform>"
    
and then run (make sure you activated the `cell-factory-design-course` environment first)

    pip install .


## Download course materials

 You can download all course notebooks at once [here](https://github.com/agilebiofoundry/cell-factory-design-course/archive/master.zip).
