
### Requirements

Three components are needed beyond MATLAB itself (R2024b or later):

1. **A local BOPTEST server** running the `multizone_residential_hydronic` testcase, deployed as a Docker container and reached over its HTTP interface at `http://127.0.0.1:80`. All building simulation happens inside it.
2. **The OSQP solver**, used through its prebuilt MATLAB interface. The interface is shipped in `osqp/` for 64-bit Windows and has to be obtained separately on any other platform.
3. **Two MATLAB toolboxes**: the Control System Toolbox, for the observer design of Chapters 3 and 4, and the Statistics and Machine Learning Toolbox, whose `bayesopt` routine implements the Bayesian optimisation of Chapter 6.


### Code layout

| Directory | Contents |
|---|---|
| `boptest/` | HTTP client wrapping the BOPTEST API, closed-loop runners of every controller, RL and BO tuning scripts, and the forecast-uncertainty and out-of-sample sweeps of Chapter 7 |
| `identification/` | Excitation and estimation code of Chapter 3 |
| `mpc/` | Matrix builder of Appendix A, the MPC step functions, the OSQP wrapper, and the helpers that compute the standard KPI report and comfort figure of every run |
| `data/` | Identification and validation datasets collected from the testcase |
| `models/` | Every artefact the pipeline produces: the identified parameters, the precomputed MPC matrices, and one `_simdata.mat` file per closed-loop run |
| `figures/` | Figures drawn by the plotting scripts, in the same form in which they appear in the thesis |
| `osqp/` | Prebuilt MATLAB interface of OSQP (version 0.6.2, 64-bit Windows) |

### Running the pipeline

From the repository root:

```matlab
startup                     % adds the project directories to the MATLAB path
ekf_sysid_2C                % grey-box identification of Chapter 3 (run once)
tracking_mpc_params_build   % constant QP matrices of Appendix A (run once)
run_baselines               % menu from which every controller of the thesis is launched
```

`models/` already contains the output of the second and third steps, together with every run reported in the thesis, so those two steps only need to be repeated in order to rebuild them. For the same reason, options 13, 18, 19 and 20 of `run_baselines` work without a BOPTEST server, since they only read `.mat` files. Every other option needs the server. Each benchmark evaluates the parameters saved by its tuner: option 5 reads what option 4 saved, 7 what 6 saved, 9 what 8 saved, and 11 what 10 saved.

### Script-to-result map

Every script saves the official BOPTEST KPIs alongside its logs.

| Script | Produces | Result in the thesis | `run_baselines` option |
|---|---|---|---|
| `identification/ekf_sysid_2C.m` | 2RC parameters | Chapter 3 results | setup step |
| `mpc/tracking_mpc_params_build.m` | QP matrices | Appendix A | setup step |
| `boptest/pi_baseline.m` | PI run | PI rows, Tables 4.3 and 7.4 | 1 |
| `boptest/mpc_baseline.m` | MPC run | MPC row, Figure 4.1 | 2 |
| `boptest/ekf_mpc.m` | Offset-free run | Offset-free row, Figure 4.2 | 3 |
| `boptest/rl_mpc_loop.m` | RL training episode | Chapter 5 training results | 4, 6 |
| `boptest/rl_mpc_benchmark.m` | Greedy RL evaluation | RL rows, Tables 5.1 and 7.4 | 5, 7 |
| `boptest/bo_mpc_tune.m` | BO search history | Figure 6.1, Table 6.1 | 8, 10 |
| `boptest/bo_mpc_benchmark.m` | BO-tuned run | Table 6.2, Figure 6.2 | 9, 11 |
| `boptest/sweep_uncertainty.m` | Degraded-forecast runs | Table 7.4 | 12 |
| `boptest/build_uncertainty_report.m` | Comparison figures | Figure 7.1 | 13 |
| `boptest/sweep_out_of_sample.m` | Unseen-week runs | Table 7.5 | 14 to 17 |
| `boptest/plot_out_of_sample.m` | Unseen-week figure | Figure 7.2 | 18 |
| `boptest/gradient_check.m` | Finite-difference check | Section 5.4 | 20 |
| `generate_controller_figures.m` | Controller figures | Chapters 4 to 6 | 19 |

Options 10 and 11 repeat the Bayesian optimisation on top of the offset-free MPC, which is the control experiment reported in Section 6.2.

### Datasets

The data-collection scripts in `boptest/` regenerate the datasets of `data/` from the testcase. `data_collection_mlprs.m` applies the multilevel excitation of Chapter 3 and writes `multilevel_identification.csv` and `multilevel_validation.csv`, while `data_collection_baseline.m` records the closed-loop operation of the built-in PI controllers and writes `baseline_identification.csv` and `baseline_validation.csv`. These scripts only need to be re-run when a new dataset is required, since the collected CSV files are kept in the repository.
