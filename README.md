# Optimal Subsampling Designs

Code used for the experiments in the paper Optimal Subsampling Designs by Henrik Imberg, Marina Axelson-Fisk, and Johan Jonasson. 

* main.R: Main function for the experiments in Section 6.

* acov.R: Evaluates the approximate covariance matrix (Eq. (7)-(9), Section 2.2).

* Phi.R: Evaluates the Phi-optimality objective function (Table 1, Section 2.3). 

* optimise_ss.R: Finds optimal sampling scheme using the iterative L-approximation algorithm (Algorithm 2, Section 3.3).

* fmt.R: Simple text formatting function used for Table 2–4 (Section 6.2–6.4).

* figures.R: Generates Figure S1 and S2 in Appendix B.

Code is included in the RScript folder.

Input data is stored in the Data folder. 

Generated tables and figures are stored in the Tables and Figures folders.

Required packages: cowplot, expm, sampling, tidyverse, and xtable.