simplesims
======================

Simulation (and analysis) code for testing brain functional connectivity methods

When using these simulations in publications please cite:
Cole, Yang, Murray, Repovs, Anticevic, "Functional connectivity change as shared signal dynamics", Journal of Neuroscience Methods, http://dx.doi.org/10.1016/j.jneumeth.2015.11.011

Contact: Michael W. Cole, mwcole@mwcole.net, http://www.colelab.org


The code can be run in R, using the following commands:

source('CorrelationVsCovarianceSim.R')

source('CorrelationVsCovarianceSim_WithOtherMeasures.R')


Before running the code, be sure to have R installed (e.g., RStudio), and that you have the following R packages installed:
ggplot2, sapa, plyr, entropy, grid


An alternate version of the simulation framework for MATLAB is also included:
CorrelationVsCovarianceSim_matlab.m

