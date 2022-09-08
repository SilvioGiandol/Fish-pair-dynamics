fivesimulations.m:matlab script to make simulations over the 5D system
ODEfive.m: matlab function associated to fivesimulations.m

tandemheatmaps.m: generate heatmaps of stable configuration for tandem
staggeredheatmaps.m: generate heatmaps of stable configuration for staggered
They need the functions:
SOLVEeq.m: equations to solve
equationmap.m: check if all the equations are satisfied
jacob.m: check the eigenvalues of the jacobian + displays eigenvalues and eigenvectors + finds frequencies

crossstreamequilibria: plot cross-stram equilibria graphs for qualitative dynamics, needs the same functions as tandemheatmaps.m