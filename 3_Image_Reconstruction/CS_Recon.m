%CS reconstruction from a undersampled MRI data. 


addpath('utilities');
addpath('solver');

%% Import Bruker MRI raw file.
npe0=128;
pn='C:\Users\Ming\Documents\MATLAB\Compressed Sensing\CSrec\Data\CS 04242013 cineheart\Fully sampled';
[FB mask m n] = readbruker(npe0,pn);

F = fftshift(fft2(FB(:,:,1),m,n));%test the first frame if cine
imF = mat2gray(abs(F)); %Magnitude of original image for display purpose

%% Modify the random mask into 'picks' as input of RecPF solver. 
pn='C:\Users\Ming\Documents\MATLAB\Compressed Sensing\CSrec\Data\CS 04242013 cineheart\CS cine';
[FB2 mask m2 n2] = readbruker(npe0,pn);

mask0 = mask;
mask = fftshift(mask0); %rearrange the mask origins as generate input of solver.
ky=n;

%Convert mask to index of sampling position 'picks' used in RecPF solver
a = [1:ky]';
maskind = mask.*a; 
maskind(maskind==0)=[];
npe = size(maskind,1);

for a=2:ky
    maskind(:,a)=maskind(:,1)+ky*(a-1);
end

picks = reshape(maskind, npe*ky, 1);%picks is the input of RecPF solver
%picks has the origin at the corner and it stores the index of sampled
%kspace point. 

%% Create the undersampled simulation data with the generated mask

[m,n]=size(F);
FB1 = ifft2(F)/sqrt(m*n);%convert original image F into kspace. 
%FB1=FB(:,:,1);
B = FB1(picks);%simulate the undersampled kspace by only picks the points indexed by picks.

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
[U,Out_RecPF] = RecPF(m,n,aTV,aL1,pick,B,2,opts,WT,W,1);
toc
disp('done!');

U2 = mat2gray(abs(U));

%% plot image

figure(1);
subplot(131); imshow(F,[]); title('Original');
subplot(132); 
bar(mask0);
%hold on, plot(pdf/pdf(ky/2+1),'r--')%bar(pick(:,1)); 
title('Undersampling pattern');
subplot(133); imshow(U2,[]); title(sprintf('Recon. RelErr=%4.2f%% SNR=%4.1f',norm(U2(:)-F(:))*100/norm(F(:)),snr(U2,F)));
