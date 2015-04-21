% This is a demonstration of a real randomly undersampled 3DFT acquisition
% It includes normalization, phase correction and separable reconstruction of each slice
% 
% In each step, the previous slice is used as an initial image, since they are very close 
% in content.
%
% The full recon takes ~ 430 seconds on an AMD64 3700+ processor with 2.5Gb of RAM
%
% The data was obtained with collaboration with Tolga Cukur.
% The original data was 324x180x180. Here only 64 slices in the middle of the calf
% are given because of the large file size
%
% The data was undersampled in the acquisition by a factor of 2
% the zero filling with density compensation, not surprizingly, result in a good MIP, but has
% significantly more "noise" the slices than the CS recon.
%

load kspaceDS.10.32.mat


% take ifft in the fully sampled dimension -- ALREADY DONE IN DATASET
im_zfwdc = zeros(size(data));
for n=1:size(data,3)
	im_zfwdc(:,:,n) = ifft2c(data(:,:,n)./pdf);
end

% scale data such that the maximum image pixel in zf-w/dc is around 1
% this way, we can use similar lambda for different problems
data = data/max(abs(im_zfwdc(:)));
im_zfwdc = im_zfwdc/max(abs(im_zfwdc(:)));