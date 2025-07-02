# Global Sensitivity Analysis using LHS-PRCC

## Overview
This R script performs a global sensitivity analysis on a 16-state epidemiological model using Latin Hypercube Sampling (LHS) and Partial Rank Correlation Coefficients (PRCC). It evaluates the impact of parameters on key outcomes, such as cumulative resistant infections and deaths, over a 10-year simulation period.

## Dependencies
The script requires the following R packages:
- `deSolve`: Solves ordinary differential equations (ODEs).
- `lhs`: Generates Latin Hypercube Samples.
- `sensitivity`: Computes PRCC for sensitivity analysis.
- `dplyr`: Handles data manipulation.
- `ggplot2`: Creates visualizations (e.g., PRCC heatmap).
- `doParallel`: Enables parallel processing.
- `purrr`: Supports functional programming.

Install these packages using:
```R
install.packages(c("deSolve", "lhs", "sensitivity", "dplyr", "ggplot2", "doParallel", "purrr"))
```

## Script Structure
1. **Load Libraries**: Imports required R packages.
2. **Helper Function**: Defines `get.beta.par` to fit Beta distribution parameters from quantiles.
3. **Parallel Processing**: Configures parallel computing to optimize simulation speed.
4. **Model Setup**: Defines initial conditions, simulation time (10 years), and the ODE model (`infection_model_ode`) with constrained flows for numerical stability.
5. **Parameter Sampling**: Uses LHS to generate 1000 parameter sets from Beta and Uniform distributions based on provided quantiles.
6. **Simulations**: Runs the ODE model in parallel for each parameter set.
7. **PRCC Calculation**: Computes PRCC for outcomes (e.g., cumulative hospital/community resistant infections, deaths).
8. **Visualization**: Generates a heatmap of PRCC results to visualize parameter impacts.

## Key Features
- **Model**: 16-state epidemiological model with hospital and community dynamics.
- **Sampling**: LHS for efficient parameter space exploration (1000 samples).
- **Outcomes**: Focuses on cumulative resistant infections (`C_RIh`, `C_RIc`) and deaths (`Drh`, `Drc`).
- **Parallelization**: Uses `doParallel` for faster computation.
- **Visualization**: PRCC heatmap with `ggplot2` for clear sensitivity analysis results.

## Usage
1. Ensure all dependencies are installed.
2. Run the script in an R environment (e.g., RStudio).
3. The script will:
   - Generate parameter samples using LHS.
   - Run simulations in parallel.
   - Compute PRCC for key outcomes.
   - Display a PRCC heatmap (optionally saved as `prcc_heatmap.png` if uncommented).
4. Check console output for progress and simulation success rate.

## Notes
- **Reproducibility**: Set seed (`set.seed(123)`) ensures consistent results.
- **Simulation Count**: Default is 1000 samples; adjust `n_samples` for precision vs. speed.
- **Output**: PRCC heatmap shows parameter influence on outcomes; values range from -1 to 1.
- **Error Handling**: Failed simulations are filtered out; check console for issues.
- **Customization**: Modify `param_data` for different parameter ranges or distributions.

## Output
- **Console**: Logs progress (e.g., number of successful simulations).
- **Plot**: PRCC heatmap visualizing parameter sensitivities.
- **Optional Save**: Uncomment `ggsave` to save the heatmap as a PNG file.

## Contact
For issues or questions, contact Yewei Xie at yewei.xie@u.duke.nus.edu.


