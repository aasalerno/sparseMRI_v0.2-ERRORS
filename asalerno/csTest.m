close all;
clear all;
file = {'/projects/muisjes/asalerno/CS/data/RealImgRaw.10.28.mnc' ...
    '/projects/muisjes/asalerno/CS/data/ImagImgRaw.10.28.mnc'};
sty = 'circ'; % Fully sampled region style
sampFac = 0.5; % Undersampling factor
sl = 250; % slice that we want to get
loc = 1;
[im,fil] = testMap(file,sty,sl,loc);
data = fft2(im);

N = size(im);		% image Size
DN = N;         	% data Size
pctg = sampFac;  	% undersampling factor
P = 5;              % Variable density polymonial degree
TVWeight = 0.1; 	% Weight for TV penalty
xfmWeight = 0.1;	% Weight for Transform L1 penalty
Itnlim = 30;         % Number of iterations

% initialize Parameters for reconstruction
phmask = zpad(hamming(6)*hamming(6)',N(1),N(2)); %mask to grab center frequency
phmask = phmask/max(phmask(:));			 %for low-order phase estimation and correction
ph = exp(1i*angle((ifft2c(data.*phmask))));



% generate variable density random sampling
pdf = genPDF(DN,P,pctg,2,0.1,0);	% generates the sampling PDF - L2
% In this case, there is a clear dropoff in the density of the sampling as
% one gets away from the centre - optially it appears as if there is an
% isotropic fall off
k = genSampling(pdf,10,60);		% generates a sampling pattern
k = k | fil;
% Generate the required operators
FT = p2DFT(k,N,ph,2);
XFM = Wavelet('Daubechies',6,4);	% Wavelet - gives a 1x1 matrix, of value 1
%XFM = TIDCT(8,4);			% DCT
%XFM = 1;				% Identity transform


params = init;
params.FT = FT;
params.XFM = XFM;
params.TV = TVOP;
params.data = data;
params.TVWeight =TVWeight;     % TV penalty
params.xfmWeight = xfmWeight;  % L1 wavelet penalty
params.Itnlim = Itnlim;


im_dc = fftshift(FT'*(data./pdf));	% init with zf-w/dc (zero-fill with density compensation)
res = XFM*im_dc; % This is the undersampled image


x0 = res; % Undersampled

tic
for n=1:Itnlim
    res = fnlCg(res,params);
    im_res = XFM'*res;
    if n==Itnlim; imshow(abs(im_res),[]); end
end
toc

h = figure;
subplot(121);
imshow(flipud(abs(im)'),[]);
subplot(122);
imshow(flipud(abs(im_res)'),[]);

st = input('Save?');
if strcmp('y',st(1))
    saveas(h,['/micehome/asalerno/Dropbox/CSRecon' num2str(sl) '.jpg'])
end
% subplot(133);
% imshow(flipud(-(abs(im_res)-abs(im_dc))'),[]);
% colorbar
% suptitle('Original undersampled                        CS Recon')
% mask_lr = genLRSampling_pctg(DN,pctg,1,0); % Creates a mask just to view a low res image
% im_lr = ifft2c(zpad(data.*mask_lr,N(1),N(2)));
% 
% im_full = ifft2c(zpad(fft2c(data),N(1),N(2))); % Full image
% figure, imshow(abs(cat(2,im_full,im_lr,im_dc,im_res)),[]);
% title('original             low-res              zf-w/dc              TV');
% 
% figure, plot(1:N(1), abs(im_full(end/2,:)),1:N(1), abs(im_lr(end/2,:)), 1:N(2), abs(im_dc(end/2,:)), 1:N(2), abs(im_res(end/2,:)),'LineWidth',2);
% legend('original', 'LR', 'zf-w/dc', 'TV');
