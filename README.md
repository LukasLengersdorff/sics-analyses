## Sampling in Constrained Space: Code & Analyses

As reported in *Sampling in Constrained Space: Efficient Estimation of Model Evidence under Equality and Inequality Constraints* (Lengersdorff & Marsman, 2025).

### Setup

Copy/clone this repository to your computer. Then, you first need to install the necessary packages, for which you can use the script `install_packages.R`. This will install our prototype `sics` package, which is an R interface to SICS-compatible STAN models, as well as a few additional packages needed for our analyses and plots.
The `sics` package can be found in this repository: https://github.com/LukasLengersdorff/sics-prototype

### Overview of analyses

-   **Sampling demo:** Code to demonstrate HMC sampling in unconstrained and constrained space (Figure 2).
-   **Single Subject:** Single-subject reinforcement learning model analyses, with comparison between SICS and the encompassing prior methods (Figure 3 & 4 and numerical results).
-   **Multilevel:** Constrained multilevel reinforcement learning models (Figure 5 and numerical results).
