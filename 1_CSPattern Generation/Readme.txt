Run matlab cs_pattern.m file in the folder named cs_pattern. 

In the m file, important inputs are:
trace: the size of ky of a fully sampled k-space.
ratio: compressive ratio

This m-file calls function cs_generatemask. Change the input of this function to acquire different sampling pattern. 

The outputs of this m-file  are cs_phaseencoding.mat file and ACQ_spatial_phase1.txt. 

The .mat file contains two variables: penum as the k-space index of acquiring lines and grdpwr as the gradient power of phase encoding gradient. 

The ACQ_spatial_phase1.txt is ASCII print of grdpwr values with "," as delimiter. It will be used by the ParaVision software when generating the user defined options.

