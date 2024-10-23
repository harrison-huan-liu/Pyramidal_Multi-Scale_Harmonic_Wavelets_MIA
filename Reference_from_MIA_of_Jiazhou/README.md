# Characterizing the Propagation Pathway of Neuropathological Events of Alzheimer's Disease Using Harmonic Wavelet Analysis

This repository contains the code that was used for the harmonic wavelets identification on brain network data.

All preprocessed brain network data that were used in the analysis are available in folder "Data", including clinical information file "DataTS.cvs" and brain network file in folder "AD-Data";

"main.m" is the main program to obtain the final region-adaptive common harmonic wavelets for each brain region;
"Preprocess_network_data.m" loads and preprocesses all brain networks.
"Estimate_glob_com_harmonics.m" estimates the global common harmonic waves with the help of "Calculate_IndividualPhi.m" (iteratively update each individual harmonic waves); 
"Construct_individual_harmonic_wavelets.m" constructs the region-adaptive individual harmonic wavelets;
"Identify_glob_com_har_wavelets.m" identifies the region-adaptive common harmonic wavelets.
