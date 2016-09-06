---
layout: page
title: Setup
permalink: /setup/
---

## Prerequisites

Before installing any of the course-specific packages below, please head over to <https://biosustain.github.io/2016-09-05-software-carpentry/> and
follow the installation instructions for Python.

## Installing cameo

Please open a shell (terminal) and create a conda environment for the course using the following command.

    conda create -y --name cell-factory-course python=3.4 scipy pandas numexpr bokeh jupyter sympy lxml matplotlib

If you're on Linux or OS X or are using the git-bash shell on Windows run

    source activate cell-factory-course

to activate the environment. If you're on Windows using the default terminal (cmd) run the following command instead.

    activate cell-factory-course

If you successfully activated you're environment, your command prompt will look similar to this.

    (cell-factory-course)$

Then use `pip` to install [cameo](http://cameo.bio) and [driven](http://driven.bio).

    pip install cameo driven

Run the following command to check that the installation was successful.

    python -c "from cameo import models;print(models.bigg.e_coli_core.solve().objective_value)"

The output should be `0.8739215069684304`.

## Troubleshooting

If you're on a Mac, the test command might fail because a solver was not properly installed. Please install homebrew

## Download course materials

 You can download all course notebooks at once [here](https://github.com/biosustain/cell-factory-design-course).