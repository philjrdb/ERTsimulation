# ERTsimulation

Code used in Jean-Richard-dit-Bressel et al. (2020). "Analyzing Event-Related Transients: Confidence Intervals, Permutation Tests, and Consecutive Thresholds", Front. Mol. Neurosci. https://doi.org/10.3389/fnmol.2020.00014


*ERT_simulation.m* generates 2 artificial peri-event waveform populations (one with a temporally-aligned transient [signal/ERT popn] and one without [control/null popn]). It then runs a Monte Carlo simulation of analysis methods (t confidence intervals, bootstrapped confidence intervals and permutation test), with and without a consecutive threshold. Outputs results in command window, diary log & figures.

Parameters can be edited at top of script as needed. Uses in-built *lowpass* function (needs MATLAB 2018a or newer).

Various analysis/report functions:
*boot_CI.m* , *boot_diffCI.m* , *consec_idx.m* , *ERT_sim_report.m* , *permTest_array.m*, *sem.m*, *CIadjust.m*, *errorplot3.m*

*Phil_ERTrealexample.m* === script that applies same analyses to real data (ERT_testdata.mat), outputs figures


***********************************************************************
EDITS/ADDITIONS:
[2020-04-17] 
Added missing functions *sem.m* *CIadjust.m* *errorplot3.m*
Removed unnecessary data-import code from *Phil_ERTrealexample.m* , added code to label significance bars

***********************************************************************
EDITS/ADDITIONS:
[2020-11-02] 
Added missing function *col_rep.m* (colour repository for graphs)
