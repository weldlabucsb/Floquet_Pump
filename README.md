# Floquet_Pump
Theory and analysis for Lithium Floquet Pump Experiment

This repository is for theory simulations and data analysis code related to the Floquet pumping, the use of lattice amplitude modulation
for pumping atoms up and down a confining harmonic potential.

**psbo_semi**: In this folder, there are semiclassical simulations of the Floquet pumping scheme. Grab the whole folder to run these
simulations. Keep the bands.fig file in this folder as it is used to compute band dispersions. Run each script while in the folder of that version. The current version control is:
1. Single particle position evolution, good and fast for first moment analysis.
2. Classical ensemble first and second moment evolution, there are some issues with the momentum calculation.
3. Same as v2 with parallel computing and standardized scale for integration time step.
4. Added probability distribution evolution (OD) plot by binning points in position space and interpolating. Fixed error in momentum moments calculation by projecting to 1st Brillouin zone prior to computing moment. Now produces a data.mat file for comparison with experiment.
5. Added correction for half Bloch oscillation that occurs prior to the start of the experiment by computing a new initial force assuming a Wannier Stark length of motion from the inputted initial conditions. Fixed error in momentum distribution initialization that was causing uncertainty to be less than hbar/2.
6. Now uses a binned probability distribution to calculate the position moments; bin size is based on diffraction limit of imaging system.

**psbo_quantum**: In this folder is a quantum simulation of the Floquet pumping scheme. Grab the whole folder and run the script Bloch oscillate. It time evolves an initial Bloch wave with a Gaussian width in position space. Parameters are set in two sections, experimental (drive frequency, Bloch frequency etc.) and simulation parameters (mesh size, time step). It can do either Runge-Kutta 4 evolution or midpoint Crank Nicholson, though the CK is broken I think. This is set in the time evolve section, either calling makeRK4 or makeCK. Though there are normally hundreds of thousands to millions of time evolve steps, the output data is sampled only 1000 times over the interval (set by inc).

**theory_exp_compare** In this folder, there is analysis script theory_exp_compare_vN.m to compare the experimental data to the semiclassical simulations of psbo_semi. To use, create a folder inside that contains the corresponding atomdata.mat from the experiment and data.mat from psbo_semi. Then run the script. Current version control is:
1. Computes first and second moments of experiment by summing a weighted OD. Generates a plot overlaying the theory and experiment for the first and second moments.
