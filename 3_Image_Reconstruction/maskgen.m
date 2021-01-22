function [mask pdf] = maskgen(ky,t,ratio,samplingtype,var,symmetry)
% maskgen creates random under-sampled phase encoding scheme for compressive sensing MRI scan. It writes the gradient power array into ACQ-phase-create time.txt file.
%
% [mask grdpwr] = maskgen(ky,t,ratio,samplingtype,symmetry) generates t frames of phase encoding scheme undersampled to a ratio of original ky-size scheme. The undersampling pattern follows random distribution of samplingtype with center symmetry (symmetry=1) or asymmetry (symmetry =0) 
%
% mask is a ky*t output matrix contains 0 and 1. It has value 0 for not sampling corresponding phase encoding step. 
%
% grdpwr is a cell array of t vectors. Each element of grdpwr stores the phase encoding gradient input as ACQ_Spatial_Phase1 in Bruker paravision5.1 corresponding to mask(:, t). 
%
% ky is a integer, indicating the number of phase encoding steps of fully sampled kspace. 
%
% t is  a integer, indicating the number of frames of video/kspaces.
%
% ratio is a double between 0 and 1, indicating  the compressive ratio. 
%
% samplingtype is the random distribution function. 
%                          'gaussian'          normal distribution
%
% var is the parameter other than mean used in the pdf of input sampling type.
%
% symmetry  indicates whether the phase encoding scheme is center asymmetric. 
%                          0             asymetry
%                          1             symmetry
%
%   $Author: Ming Yang my5f2@mail.missouri.edu
%   $version: 1.1 $  $Date: 08-Jul-2013 16:02:07 $
%
% new case added to samplingtype: Cauchy (Lorentz) distribution
%
%   $modified by: Elliot Tan et7@princeton.edu
%   $version: 1.2 $  $Date: 09-Jul-2013 14:51:20 $
% commented by Rosa Zheng 2013_0709
%
% new case added to samplingtype: uniform distribution
% added other possibilities in Cauchy and uniform case
%
%   $modified by: Elliot Tan et7@princeton.edu
%   $version: 1.3 $  $Date: 10-Jul-2013 12:49:42 $
% new case added to samplingtype: uniform distribution
% added other possibilities in Cauchy and uniform case
%
%   $modified by: Elliot Tan et7@princeton.edu
%   $version: 1.3 $  $Date: 10-Jul-2013 12:49:42 $
%
% Employed a cell code structure and changed the output Files to a certain
% folder. 
% Fixed the bug when symmetry=1. 
%   $modified by: Ming Yang my5f2@mail.missouri.edu
%   $version: 1.4 $  $Date: 12-Jul-2013 17:03:07 $

%% Generate PDF array according to input options
% 
% %for debug
% clear all; close all; clc;
% ky =256;t=3e3;ratio=0.2;var =21;samplingtype='G';symmetry=0;
% % comment out the above after debugging

 mean = ky/2; % moved to beginning, by Rosa Zheng
 a = 1:ky; % use vector instead of for-loop
 PI = 3.141593;
switch samplingtype
    case {'Gaussian','gaussian', 'G','g'} % by Rosa Zheng
         const1 = sqrt(2*pi);
         const2 = 2*sqrt(2);
         pdf(a) = exp(-(a-ky/2).^2/(2*var^2));%true Gaussian pdf =exp(-(a-ky/2).^2/(2*sigma^2))/(const1*sigma);        
%calculate the scale for discrete R.V. X, so that the sum(P(X<=pdf))/(ky*scale)=ratio.
         scale = const1*var*erf(ky/(const2*var))/(ky*ratio);
        
    case {'Cauchy', 'cauchy', 'C', 'c'} % pls add other possibilities as in the Gassian case -Rosa Zheng
% Specifically for Cauchy dist.:
% CDF = (1/pi) arctan[(x - x0)/gamma] + .5
% PDF = 1 / [pi*gamma(1 + ((x - x0) / gamma)^2)] 
		
%Generating scaled PDF for Cauchy dist. using solved gamma variable	
         pdf(a) = 1 ./ ( (1 + ((a - mean) / var) .^2)); %true  pdf= 1 ./ (gamma*pi* (1 + ((a - mean) / gamma) .^2));  
%calculate the scale for discrete R.V. X, so that the sum(P(X<=pdf))/(ky*scale)=ratio.
         scale = (var * atan(ky / var)) / (ky*ratio);
         
    case {'vonMises', 'vonmises', 'Vonmises', 'V', 'v'}
% Generating the pdf of von Mises
% Original pdf function of von Mises: pdf(a) = (exp(var * cos(x - mean))) / (2 * PI * besseli(0, var)) 
% where var = k, which is a measure of concentration, or 1/variance   

    b = linspace(0, 2 * PI, ky); % Vector needed since support of von Mises is [0, 2*PI], mean = PI
	pdf(a) = (exp(var * cos(b - PI))) / (2 * PI * besseli(0, var)); % Setting probabilities of original support [0, 2*PI] into new support of [1, ky]

% calculate the scale for discrete R.V. X, so that the sum(P(X<=pdf))/(ky*scale)=ratio
% Need to satisfy equation ratio = (CDF(ky) - CDF(0))/(scale * ky)
% CDF(ky) - CDF(0) = integral of pdf(x) from 0 to ky
    
    % Riemann sum method
    dx = .001;
    range = 0:dx:ky;
    scale = sum(((exp(var * cos(range - PI))) / (2 * PI * besseli(0, var))) * dx) / (ky * ratio);
    
    case {'Uniform', 'uniform', 'U', 'u'}
         pdf(a)=ratio;
	
    otherwise 
        sprintf('PDF is not supported');
end


%% Implement Monte-Carlo Method to generate discrete R.V. array 'mask' following PDF.
mask = zeros(ky,t);

for indt = 1:t       
    switch symmetry
    case 0
        for a=1:ky
            if (pdf(a)>=rand*scale) %/(const1*sigma)) % this is inverse transform method
                % comments by Rosa Zheng
                % reference: Simulation, chapter 4, by Sheldon M. Ross
                % scale rand to satisfy the sampling ratio. 
                mask(a,indt)=1;
            end
        end
    case 1
        for a=1:0.5*ky
            if (pdf(a)>=rand*scale)
                mask(a,indt)=1;
            end
        end
        mask(ky/2+1:ky,indt) = mask(1:ky/2,indt)'*fliplr(eye(ky/2)); % why so cmplicated?
        
    otherwise 
        sprintf('input of symmetry option must be 1 or 0 ')
    end
    mask(ky/2+1,indt)=1;% It is already ensured, knowing pdf(ky/2)=1;
end

%Calculate grdpwr from mask.
scale = double(2/ky);
grdpwr = cell(t,1);
for indt = 1:t
    grdpwr{indt} = mask(1:ky,indt).*scale.*[-ky/2+1:ky/2]';
    for a = 1:ky
        if(grdpwr{indt}(a)==0 && a~=ky/2)
            grdpwr{indt}(a)=2; %mark the unwanted elements as 2
        end
    end
    grdpwr{indt}(grdpwr{indt}==2) = []; %delete all unwanted elements.
    grdpwr{indt} = grdpwr{indt}';
end

%% save the matrice and print grdpwr in a .txt file
filename =[ 'Sampling Scheme\' sprintf('ACQ-phase-ky%d-size%d',ky,size(mask(mask==1),1))];
save([filename '.mat'],'mask','grdpwr');
FID = fopen([filename '.txt'],'w');

grdtext = '';
for indt = 1:t
    grdtext = [sprintf('%6.6f,',grdpwr{indt}) '\n'];
    fprintf(FID, grdtext);
end

fclose(FID);

% %% Plot for verification - added by Rosa Zheng
% figure, hold on
% for indt=1:min(3,t)
%     stem(mask(:,indt)/indt)
% end    
% 
% hist_ky=sum(mask,2)/t;
% figure
% plot(hist_ky)
% hold on, plot(pdf/pdf(ky/2+1),'r--')
%             