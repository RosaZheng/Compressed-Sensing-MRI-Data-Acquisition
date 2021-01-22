clear all; close all; clc;

%
% this program generates the random undersampling  patterns for MRI scanner
% lines in ky are randomly selected with Gaussian pdf to reduce sampling
% and achieve a given ratio

traces = 128; %the size of ky of a fully sampled k-space.
samplingdensity = 1.5;
ratio = 0.25; % percentage of undersampled k-spce lines for compressive sensing

[mask npdf] = cs_generatemask([1 traces], ratio, 'uniform', 'Incoherent', 'Line', 'symmetric', samplingdensity, 100);
scal = double(2/traces);
mask(traces/2+1)=1;%ensure the center line is collected.

count=sum(mask);
for i=1:traces/2
    pestepup(i) = mask(65-i).*scal.*(i);
end

for i=1:64
    pesteplow(i) = -1.*pestepup(65-i);
end

pestep = [pesteplow,0,pestepup];

grdpwr = pestep;
grdpwr(grdpwr==0)=[];

penum=find(pestep~=0);length(penum);

save('cs_phaseencoding.mat','pestep','penum');
csvwrite('ACQ_Spatial_Phase1.txt',grdpwr);
