function gradObj = gOBJ(x,params)
% computes the gradient of the data consistency

% Here we have a few things to look at:
% x - this is the data multiplied by the wavelet, so multiplying by the
%     transpose brings it back to normal - this is built from data, but
%     divided by the pdf
% data - this is the k-space dataset - undersampled
% FT - partial Fourier operator, will convert the data to k-space
% 
% Beauty of how this is done is that the "transpose" is actualy the
% inverse, which is likely how this is coded in the @p2DFT etc. 
% WE MUST INCORPORATE THIS

	gradObj = params.XFM*(params.FT'*(params.FT*(params.XFM'*x) - params.data));

% 1) Convert sparsified image back to normal, and put it in k-space. Find
% the difference
% 2) Put back into image space, then back into sparse space
gradObj = 2*gradObj ;