# ERTsimulation

*ERT_simulation.m* generates 2 artificial peri-event waveform populations (one with a temporally-aligned transient [signal/ERT popn] and one without [control/null popn]). It then runs a Monte Carlo simulation of analysis methods (t confidence intervals, bootstrapped confidence intervals and permutation test), with and without a consecutive threshold. Outputs results in command window, diary log & figures.

Parameters can be edited at top of script as needed. Uses in-built *lowpass* function (needs MATLAB 2018a or newer).

*boot_CI.m* , *boot_diffCI.m* , *consec_idx.m* , *ERT_sim_report.m* , *permTest_array.m* === necessary analysis/report functions 

*Phil_ERTrealexample.m* === script that applies same analyses to real data (ERT_testdata.mat), outputs figures
