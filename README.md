# Floquet_Pump
Theory and analysis for Lithium Floquet Pump Experiment

This repository is for theory simulations and data analysis code related to the Floquet pumping, the use of lattice amplitude modulation
for pumping atoms up and down a confining harmonic potential.

**psbo_semi**: In this folder, there are semiclassical simulations of the Floquet pumping scheme. Grab the whole folder to run these
simulations. Keep the bands.fig file in this folder as it is used to compute band dispersions. Run each script while in the folder of that version. The current version control is:
1. Single particle position evolution, good and fast for first moment analysis.
2. Classical ensemble first and second moment evolution, there are some issues with the momentum calculation.
3. Same as v2 with parallel computing and standardized scale for integration time step.
4. Classical ensemble with first moment, second moment and probability distribution evolution. Fixed momentum moments calculation.
