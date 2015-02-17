function [res, obj, RMS] = objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx, x,dx,t, params)
%calculated the objective function
% FTXFMtx is the FT of x (sparsified dataset, back to x-space, then FT)
% FTXFMtdx is the FT of dx (the difference between x and the non sparsified
%          set)
p = params.pNorm;

obj = FTXFMtx + t*FTXFMtdx - params.data; % Checking the difference between the sparsified data, 
                                          % plus a constant*the gradient
                                          % data, minus the original
                                          % dataset. In theory, this should
                                          % just be the contributions from
                                          % every other system in wGradient
obj = obj(:)'*obj(:); % Converts this to a single value, the sum of each value squared.

if params.TVWeight
    w = DXFMtx(:) + t*DXFMtdx(:); % w is just what the total variation will be
    TV = (w.*conj(w)+params.l1Smooth).^(p/2); % The norm is the L1 or L2
else
    TV = 0;
end

if params.xfmWeight
   w = x(:) + t*dx(:); 
   XFM = (w.*conj(w)+params.l1Smooth).^(p/2);
else
    XFM=0;
end



TV = sum(TV.*params.TVWeight(:));
XFM = sum(XFM.*params.xfmWeight(:));
RMS = sqrt(obj/sum(abs(params.data(:))>0));

res = obj + (TV) + (XFM) ;