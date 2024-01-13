# DMN-taskswitch-fMRI-code

## Code section
This section contains code relevant to the data analysis of fMRI experiment 'External task switches activate default mode regions without enhanced processing of the surrounding scene' which is submited for review. The pre-reg can be found here: osf.io/e3ntu

-/stimuli contains the code to generate the experiment. -/stimuli/pilot_practice.m contains script used for teaching and training the participants before the scan. -/stimuli/mri_pilot.m contains the script that generates the experiment as participants see it on the screen inside the scanner. 

-/behav_analysis is the Matlab script used to analyse the behavioural data, including stats on the reaction time for regular and context detection trials, and performance on the memory task post-scan. This script also generates the figures used in the scientific poster presented at OHBM 2023.
The folder -/figure_scripts generate all the figures in the paper.

-/fmri_analysis contains the Matlab master scripts (SPM required) that preprocesses the fMRI data for all subjects (including spatial realignment of functional volumes, slice-time correction, co-registration to the T1-weighted structural, and normalization to the Montreal Neurological Institute (MNI) template brain), and the script for univariate and multivariate analysis. 
