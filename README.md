# LinCombSNBmax
Linearly Combining Biomarkers by Direct Maximizing of Net Benefit

This repository has the in progress R code for the direct maximization of net benefit method proposed by Mishra et al (2019), which will be available through the University of Washington- Departmet of Biostatistics Department online dissertation repository shortly.   

The purpose of this code is to estimate linear combinations of risk markers (to produce a risk score) by directly maximizing the quantity standardized net benefit. This is a measure of clinical utility for risk score that predicts a binary outcome. 

Background.pdf gives the problem statement and what we are trying to estimate
maxSNB.R: Contains main R codes for optimization problem
sim.R: Contains R codes for simulating data 

If maxSNB.R is run in full followed by sim.R the first set of simulations should be produced.  

