************************************************************************
          RecPF v2.2: a solver for CS Reconstruction from Partial Fourier data
************************************************************************

Copyright (C) 2009-2010 Junfeng Yang, Wotao Yin, and Yin Zhang

1). Get Started
===================

   Run "demo_RecPF" or "demo_RecPF_noisy" to see how to use "RecPF" to reconstruct MR images from partial Fourier data.

   If an error occurs, please run "compile_first" to compile the C files under utilities


2). Introduction
====================
   
   RecPF refers to "reconstruct MR image from partial Fourier data".
Mathematically, the input argument B consists of measurements

                   B = F_p*ubar + omega,                              (1)

where ubar is an original m-by-n MR image, F_p*ubar is a subsample of
the Fourier transform F applied to ubar, and omega is random noise.

RecPF reconstructs ubar as the minimizer of the following TVL1-L2 model:

          min_u aTV*TV(u) + aL1*||Psi'u||_1 + .5*||F_p*u - B||^2,     (2)

where aTV, aL1 > 0 are regularization parameters and Psi is a sparsifying 
basis (Psi'u is either sparse or compressible). One can set aL1 = 0.

IMPORTANT notices:

1. The underlying Fourier transform should be normalized. However, MATLAB's
Fourier transforms (e.g., fft2) are not. To normalize fft2, use fft2( )/sqrt(m*n)
for m-by-n images.

2. TV(u) is the total variation of u. A standard grid with unit-length edges
is used. Therefore, the **isotropic** version of TV(u) is defined as

TV(u) = sum_{all pixels (i,j)} ( |u(i,j)-u(i+1,j)|^2 + |u(i,j)-u(i,j+1)|^2 )^(1/2).

The 4-neighbor **anisotropic** version of TV(u) is defined as

TV(u) = sum_{all pixels (i,j)} |u(i,j)-u(i+1,j)| + |u(i,j)-u(i,j+1)|.

Using non-homogeneous grids requires adding weights to the finite differences.

3. The periodic boundary conditions are applied, i.e., images are duplicated
outside their original domain (typically, rectangles).

4. The main program RecPF.m can optionally normalize the input arguments
   U, aTV, and aL1 in (2). The switch is opts.normalize.

If opts.normalize is set to true, then the following are applied

B <-- B / Greyscale_Range_of_U
aTV <-- aTV * number_of_measurements / sqrt(m*n)
aL1 <-- aL1 * number_of_measurements / sqrt(m*n).

They make parameters rather independent of the image size, pixel intensity range,
and number of CS measurements.


3). Details about "RecPF.m" 
=============================

RecPF is used as:
    
       [U,Out] = RecPF(m,n,aTV,aL1,picks,B,TVtype,opts,PsiT,Psi,URange,uOrg);

where
     aTV, aL1  --- parameters in model (2);
           IMPORTANT: see Lines 60, 66, 67 for data/parameter normalization
     picks     --- a vector with its components the indices where Fourier 
                   coefficients are taken;
     B         --- partial Fourier data, a vector;
     TVtype    --- 1 or 2, corresponding to anisotropic/isotropic TV in (2), recommend 2;
     opts      --- contains parameters for algorithm
                 * opts.maxItr: maximal total iteration number
                 * opts.gamma: 1.0 ~ 1.618
                 * opts.beta:  1 ~ 100
                 * opts.relchg_tol: stopping tolerance of norm(U-U_previous,'fro')/norm(U,'fro')
                 * opts.normalize: whether or not normalizes parameters and data
                 * opts.real_sol: whether returns solution of real type (i.e., U = real(U))
     PsiT      --- PsiT
     Psi       --- orthonormal transform Psi
     URange    --- grayscale range of image U, e.g., 1, 255, 65535
     uOrg      --- (optional) input of the true image, for info display only 

M-files in "utilities":
  
     dctPhi        --- DCT as sparsifying basis;
     identityPhi   --- Identity operator;
     MRImask       --- used to generate mask where Fourier coefficients are taken;
     snr           --- to compute Signal-to-Noise ratio;
     Wavedb1Phi    --- Haar wavelet transform;
     funcval       --- to compute function values of objective in (2).



4). Reference
====================

    For algorithmic details, such as continution on penalty parameters and 
 optimality conditions, see references:

    J. Yang, W. Yin, and Y. Zhang, "RecPF v2: a solver for CS Reconstruction from Partial Fourier data." URL: http://www.caam.rice.edu/~optimization/L1/RecPF/.

    J. Yang, Y. Zhang and W. Yin, "A fast TVL1-L2 minimization algorithm 
for signal reconstruction from partial Fourier data", Tech. Report 08-27, 
CAAM, Rice University.

6). Contact Information
=======================

RecPF is available at: http://www.caam.rice.edu/~optimization/L1/RecPF/

Comments or suggestions? Please feel free to e-mail:

Wotao Yin, CAAM, Rice University, wotao.yin@rice.edu


7).  Copyright Notice
=======================

RecPF v2 is free software; you can redistribute it and/or modify it under 
the terms of the GNU General Public License as published by the Free 
Software Foundation; either version 3 of the License, or (at your option) 
any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details at
<http://www.gnu.org/licenses/>. 