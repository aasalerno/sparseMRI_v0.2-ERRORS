file = {'/projects/muisjes/asalerno/CS/data/RealImgRaw.10.18.mnc' ...
    '/projects/muisjes/asalerno/CS/data/ImagImgRaw.10.18.mnc'};
sty = 'circ'; % Fully sampled region style
sampFac = 0.5; % Undersampling factor
loc = 100; % slice that we want to get

[data,k,pdf] = testMap(file,sty,sampFac,loc);

N = size(data);		% image Size
DN = N;         	% data Size
pctg = sampFac;  	% undersampling factor
P = 5;              % Variable density polymonial degree
TVWeight = 0.01; 	% Weight for TV penalty
xfmWeight = 0.01;	% Weight for Transform L1 penalty
Itnlim = 8;         % Number of iterations

% Generate the required operators
FT = p2DFT(k,N,1,2);
XFM = Wavelet('Daubechies',6,4);	% Wavelet - gives a 1x1 matrix, of value 1

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
res = XFM*im_dc; % This is the undersampled image


x0 = res; % Undersampled
params = param;


tic
for n=1:Itnlim
    res = fnlCg(res,param);
    im_res = XFM'*res;
    figure(100), imshow(abs(im_res),[]), drawnow
end
toc

mask_lr = genLRSampling_pctg(DN,pctg,1,0); % Creates a mask just to view a low res image
im_lr = ifft2c(zpad(im.*mask_lr,N(1),N(2)));

im_full = ifft2c(zpad(fft2c(im),N(1),N(2))); % Full image
figure, imshow(abs(cat(2,im_full,im_lr,im_dc,im_res)),[]);
title('original             low-res              zf-w/dc              TV');

figure, plot(1:N(1), abs(im_full(end/2,:)),1:N(1), abs(im_lr(end/2,:)), 1:N(2), abs(im_dc(end/2,:)), 1:N(2), abs(im_res(end/2,:)),'LineWidth',2);
legend('original', 'LR', 'zf-w/dc', 'TV');
