function [y mask nro npe0] = readbruker(npe0,pn)
%Acquire kspace and mask for a set of bruker rawfile
%
%Auther: Ming Yang@University of Missouri
%2012


%% read acqp file for essential parameters
fid = fopen([pn '\acqp'],'r','b');

acqp = fscanf(fid,'%c');
fclose(fid);

%acquire endien
ind = findstr(acqp,'$BYTORDA');
bytorda = acqp(ind+9);
%acquire data format
ind = findstr(acqp,'$GO_raw_data_format');
dataformat = ['int' acqp(ind+23:ind+24)];
%acquire the dimension of data matrix 
ind = findstr(acqp,'$ACQ_size');
m = ind+16; n=1; buff = '';
while acqp(m)~= '#'
    buff(n)=acqp(m);
    m=m+1;
    n=n+1;
end
buff2 = textscan(buff,'%f','delimiter');

nro = buff2{1}(1)/2;
npe = buff2{1}(2);

%% acquire compressive mask
ind = findstr(acqp,'$ACQ_spatial_phase');
m = ind; n=1;
while acqp(m) ~= '-'
    m = m+1;
end
buff='';
while acqp(m) ~= '#'
    buff(n)=acqp(m);
    m=m+1;
    n=n+1;
end
buff2 = textscan(buff,'%f','delimiter');

%npe0=128;
maskind = round(buff2{1}*npe0/2+npe0/2+1);
mask = zeros(npe0,1);
mask(maskind)=1;

% [npe,a] = size(maskind);
% %npefull = 256;
% for a=2:npe0
%     maskind(:,a)=maskind(:,1)+npe0*(a-1);
% end
% 
% maskind = reshape(maskind, npe*nro, 1);

%% Import fid file as kspace data. 
fid = fopen([pn,'\fid'],'r','b');
data =fread(fid,dataformat,0,bytorda);
fclose(fid);
img_num =size(data,1)/(2*npe*nro); %Acquire the number of frames of cine data.

%arrange the data into kspace cube (kx*ky*t)i.e.(readout*phase encoding*frame)
x = reshape(complex(data(1:2:end),data(2:2:end)),nro,img_num*npe);
y = zeros(nro,npe,img_num);
for a = 1:npe
    for b = 1:img_num
        y(:,a,b) = x(:,(a-1)*img_num+b);
    end
end

%% Reconstruct image from raw data using fourier transform 
% for a = 1:img_num
% img_mag(:,:,a) = abs(fftshift(fft2(y(:,:,a),nro,npe0)));
% end
% 
% %display img
% for a = 1:img_num
% I = mat2gray(img_mag(:,:,a));
% subplot(4,4,a);
% imshow(I);
% end