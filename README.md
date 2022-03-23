# Optimal Public Expenditure with Inefficient Unemployment: Code and Data

This repository contains the code and data accompanying the paper "Optimal Public Expenditure with Inefficient Unemployment", written by [Pascal Michaillat](https://pascalmichaillat.org) and [Emmanuel Saez](https://eml.berkeley.edu/~saez/), and published in the [Review of Economic Studies](https://doi.org/10.1093/restud/rdy030) in May 2019. 

## Paper webpage

The paper and its online appendix are available at https://pascalmichaillat.org/6/.

## Figure 3

Figure 3 is produced by the MATLAB script `figure3.m`.

+ The script first calibrates the sufficient-statistic formulas (23) and (24) to describe
the onset of the Great Recession in the United States. The calibration of the two formulas is described in section 4.
+ The script then uses formula (23) to compute optimal stimulus spending and formula (24) to compute the unemployment rate reached under optimal stimulus. The formulas are used under a range of unemployment multipliers and a range of elasticities of substitution between public and private consumption.
+ The script then produces the two panels of figure 3: `figure3A.pdf`, `figure3B.pdf`.

## Figure 4

Figure 4 is produced by the MATLAB script `figure4.m`.

+ The script first calibrates the matching model with land to US data. The model is
described in sections 2.2 and 2.4 and in online appendix A. The calibration is described in online appendix A.
+ The script then computes collections of steady-state equilibria, parameterized by
different levels of aggregate demand; these collections represent the different stages of
the business cycle. The simulation procedure is described in section 5.
+ The script compares three public-expenditure policies: G/Y is constant at 16.5%, its average value in the United States for 1990-2014; G/Y is given by sufficient-statistic formula (23) (a first-order approximation to the optimal policy); and G/Y is at its optimal level, where it satisfies equation (18).
+ Last the script produces the four panels of figure 4: `figure4A.pdf`, `figure4B.pdf`, `figure4C.pdf`, `figure4D.pdf`.

## Software

The results were obtained using MATLAB R2017a on macOS High Sierra.

## License

This repository is licensed under the [MIT License](LICENSE.md).