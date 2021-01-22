function [U] = CS_simulation(ky,t,ratio,samplingtype,var,symmetry,F)
%CS reconstruction of simulated undersampled MRI data
%  
%   $Author: Ming Yang my5f2@mail.missouri.edu
%   $version: 1.0 $  $Date: 12-Jul-2013 15:36:07 $
%
%Contains code from CSrec package by J. Yang, W. Yin, and Y. Zhang, 
%"RecPF v2: a solver for CS Reconstruction from Partial Fourier data." 
%URL: http://www.caam.rice.edu/~optimization/L1/RecPF/.


addpath('utilities');
addpath('solver');



%% Generate a random mask using input options. 


% % input options for debugging
% ky=n;t=1;ratio=0.3;var=50;samplingtype='Gaussian';symmetry=1;
% %comment above lines after debugging


[mask pdf] = maskgen(ky,t,ratio,samplingtype,var,symmetry);
mask0 = mask;
mask = fftshift(mask);
%Convert mask to index of sampling position 'picks' used in RecPF solver
a = [1:ky]';
maskind = mask.*a; 
maskind(maskind==0)=[];
npe = size(maskind,1);

for a=2:ky
    maskind(:,a)=maskind(:,1)+ky*(a-1);
end

picks = reshape(maskind, npe*ky, 1);%picks is the input of RecPF solver

%% Create the undersampled simulation data with the generated mask

[m,n]=size(F);
% F = abs(fftshift(fft2(y(:,:,1),m,n)));
% F = mat2gray(F)
FB = fft2(F)/sqrt(m*n);
B = FB(picks);

%   F is the original full sampled image data
%   [m,n] = size(F);
%   Let FB = fft2(F)/sqrt(m*n) be the 2D Fourier transform of F
%   then, B = FB(picks)+noise.
%   WT = []; W = [];
%   U0 is the back-project solution from B.


%% Input Reconstruction options.
aTV = 1e-4;   % 1e-4 is good for noisy (sigma = .01) under-sampled phantom

%---------------------
aL1 = 0;       % penalty parameter of |PsiT*U|_1

W = @(x)Wavedb1Phi(x,1);
WT = @(x)Wavedb1Phi(x,0);

% For wavelets, download Rice Wavelet Tools version 2.4, and unzip at ./rwt
% addpath('rwt');
% wav = daubcqf(2);
% W = @(x) midwt(x,wav);
% WT = @(x) mdwt(x,wav);
%---------------------

opts = [];
opts.maxItr = 100; 
opts.gamma = 1.0;
opts.beta = 10;
opts.relchg_tol = 5e-4;
opts.real_sol = true;    % return solution of real type (i.e., U = real(U))

opts.normalize = true;  % New parameter/data normalization was added to make 
                        % parameters rather independent of the image size, 
                        % pixel intensity range, and number of CS measurements.

%% Call L1-norm solver RecPF to reconstruct simulation data. 
pick=false(m,n);pick(picks)=true; % convert "picks" to a logical array "pick"

disp('RecPF is running ...');
tic
[U,Out_RecPF] = RecPF(m,n,aTV,aL1,pick,B,2,opts,WT,W,range(F(:)),F);
toc
disp('done!');


%% plot image

hist_ky=sum(mask0,2)/t;

figure(1);
subplot(131); imshow(F,[]); title('Original');
subplot(132); 
bar(mask0)
hold on, plot(pdf/pdf(ky/2+1),'r--')%bar(pick(:,1)); title('Undersampling pattern');
subplot(133); imshow(U,[]); title(sprintf('Recon. RelErr=%4.2f%% SNR=%4.1f',norm(U(:)-F(:))*100/norm(F(:)),snr(U,F)));
