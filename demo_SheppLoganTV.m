% this is a script to demonstrate the original experiment by Candes, Romberg and Tao
%
% (c) Michael Lustig 2007

rand('twister',2000);
addpath(strcat(pwd,'/utils'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1 Recon Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = [256,256]; 		% image Size
DN = [256,256]; 	% data Size
pctg = [0.33];  	% undersampling factor
P = 5;			% Variable density polymonial degree
TVWeight = 0.01; 	% Weight for TV penalty
xfmWeight = 0.01;	% Weight for Transform L1 penalty
Itnlim = 8;		% Number of iterations

% generate variable density random sampling
pdf = genPDF(DN,P,pctg, 2 ,0.1,0);	% generates the sampling PDF - L2
% In this case, there is a clear dropoff in the density of the sampling as
% one gets away from the centre - optially it appears as if there is an
% isotropic fall off
k = genSampling(pdf,10,60);		% generates a sampling pattern

%generate image with noise
im = (phantom(N(1)))  + randn(N)*0.01 + 1i*randn(N)*0.01;

%generate Fourier sampling operator
FT = p2DFT(k, N, 1, 2); % FT Sampling Operator - note that this is NOT fully sampled
data = FT*im; % This builds the UNDERSAMPLED FT of the image - i.e. we are already working in undersampled k-space

%generate transform operator

XFM = Wavelet('Daubechies',6,4);	% Wavelet - gives a 1x1 matrix, of value 1
%XFM = TIDCT(8,4);			% DCT
%XFM = 1;				% Identity transform

% initialize Parameters for reconstruction
param = init;
param.FT = FT;
param.XFM = XFM;
param.TV = TVOP;
param.data = data;
param.TVWeight =TVWeight;     % TV penalty
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;

im_dc = FT'*(data./pdf);	% init with zf-w/dc (zero-fill with density compensation)
% im_dc - data = data.*(1 - 1/pdf) --> this has the density compensation,
% which is surprisingly important
figure(100), imshow(abs(im_dc),[]);drawnow;

res = XFM*im_dc; % This is the undersampled image


% added by AS
x0 = res; % Undersampled
params = param;

% do iterations
tic
for n=1:8
    res = fnlCg(res,param);
    im_res = XFM'*res;
    figure(100), imshow(abs(im_res),[]), drawnow
end
toc

%//////////////////
% ViewPorts
%//////////////////
% create a low-res mask
mask_lr = genLRSampling_pctg(DN,pctg,1,0); % Creates a mask just to view a low res image
im_lr = ifft2c(zpad(fft2c(im).*mask_lr,N(1),N(2)));

im_full = ifft2c(zpad(fft2c(im),N(1),N(2))); % Full image
figure, imshow(abs(cat(2,im_full,im_lr,im_dc,im_res)),[]);
title('original             low-res              zf-w/dc              TV');

figure, plot(1:N(1), abs(im_full(end/2,:)),1:N(1), abs(im_lr(end/2,:)), 1:N(2), abs(im_dc(end/2,:)), 1:N(2), abs(im_res(end/2,:)),'LineWidth',2);
legend('original', 'LR', 'zf-w/dc', 'TV');




