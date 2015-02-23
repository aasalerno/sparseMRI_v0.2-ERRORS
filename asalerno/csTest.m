file = {'/projects/muisjes/asalerno/CS/data/RealImgRaw.10.18.mnc' ...
    '/projects/muisjes/asalerno/CS/data/ImagImgRaw.10.18.mnc'};
sty = 'circ'; % Fully sampled region style
sampFac = 0.75; % Undersampling factor
sl = 100; % slice that we want to get
loc = 1;
im = testMap(file,sty,sl,loc);

N = size(data);		% image Size
DN = N;         	% data Size
pctg = sampFac;  	% undersampling factor
P = 5;              % Variable density polymonial degree
TVWeight = 0.01; 	% Weight for TV penalty
xfmWeight = 0.01;	% Weight for Transform L1 penalty
Itnlim = 8;         % Number of iterations

% generate variable density random sampling
pdf = genPDF(DN,P,pctg, 2 ,0.1,0);	% generates the sampling PDF - L2
% In this case, there is a clear dropoff in the density of the sampling as
% one gets away from the centre - optially it appears as if there is an
% isotropic fall off
k = genSampling(pdf,10,60);		% generates a sampling pattern

% Generate the required operators
FT = p2DFT(k,N,1,2);
data = FT*im;
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
res = XFM*im_dc; % This is the undersampled image


x0 = res; % Undersampled
params = param;


tic
for n=1:Itnlim
    res = fnlCg(res,param);
    im_res = XFM'*res;
    if n==Itnlim; imshow(abs(im_res),[]); end
end
toc

figure;
subplot(131);
imshow(abs(im_dc),[]);
subplot(132);
imshow(abs(im_res),[]);
subplot(133);
imshow(abs(im_res)-abs(im_dc),[]);
suptitle('Original undersampled           CS            Residual')
% mask_lr = genLRSampling_pctg(DN,pctg,1,0); % Creates a mask just to view a low res image
% im_lr = ifft2c(zpad(data.*mask_lr,N(1),N(2)));
% 
% im_full = ifft2c(zpad(fft2c(data),N(1),N(2))); % Full image
% figure, imshow(abs(cat(2,im_full,im_lr,im_dc,im_res)),[]);
% title('original             low-res              zf-w/dc              TV');
% 
% figure, plot(1:N(1), abs(im_full(end/2,:)),1:N(1), abs(im_lr(end/2,:)), 1:N(2), abs(im_dc(end/2,:)), 1:N(2), abs(im_res(end/2,:)),'LineWidth',2);
% legend('original', 'LR', 'zf-w/dc', 'TV');
