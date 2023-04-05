# Optimal Subsampling Designs

Code used for the experiments in the paper Optimal Subsampling Designs by Henrik Imberg, Marina Axelson-Fisk, and Johan Jonasson. 

Code is included in the Rscript folder.

* main.R: Main function for the experiments in Section 6.

* acov.R: Evaluate the approximate covariance matrix (Eq. (7)-(9), Section 2.2).

* Phi.R: Evaluate the Phi-optimality objective function (Table 1, Section 2.3). 

* Phiopt.R: Find Phi-optimal sampling scheme using the fixed-point method (Algorithm 2, Section 3.3).

* Lopt.R: Find L-optimal sampling scheme using Algorithm 1, Section 3.3.

* fmt.R: Simple text formatting function used for Table 2–4 (Section 6.2–6.4).

* figures.R: Generates Figure S1 and S2 in Appendix B.

Input data is stored in the Data folder. 

Generated tables and figures are stored in the Tables and Figures folders.

Required packages: benchmarkme, cowplot, expm, sampling, tidyverse, and xtable.
