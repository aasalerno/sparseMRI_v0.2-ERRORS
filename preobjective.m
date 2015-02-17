function [FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx] = preobjective(x, dx, params)

% precalculates transforms to make line search cheap
%
% Params.FT is the operator dictating what part of the FT we care about -
%           it is a partial fourier operator
%
% Params.XFM is the method of calculating the transform
%            - it is a sparse space
%            - It is created depending on how we want the sparsification
%
% x - this is the undersampled image
%
% dx - this is the gradient of the undersampled image, i.e. the difference
% between the undersampled image and the normal image


FTXFMtx = params.FT*(params.XFM'*x); % Applies the partial FT to the image that has been sparsified
FTXFMtdx = params.FT*(params.XFM'*dx); % Applies the partial FT to the gradient
          % An interesting note is that there seems to be minimal low
          % frequency noise, and it is zero as a circle!

if params.TVWeight % If the TV matters, calculate it with the TV too!
    DXFMtx = params.TV*(params.XFM'*x);
    DXFMtdx = params.TV*(params.XFM'*dx);
else
    DXFMtx = 0;
    DXFMtdx = 0;
end