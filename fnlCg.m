function x = fnlCg(x0,params)
%-----------------------------------------------------------------------
%
% res = fnlCg(x0,params)
%
% implementation of a L1 penalized non linear conjugate gradient reconstruction
%
% The function solves the following problem:
%
% given k-space measurments y, and a fourier operator F the function 
% finds the image x that minimizes:
%
% Phi(x) = ||F* W' *x - y||^2 + lambda1*|x|_1 + lambda2*TV(W'*x) 
%
%
% the optimization method used is non linear conjugate gradient with fast&cheap backtracking
% line-search.
% 
% (c) Michael Lustig 2007
%-------------------------------------------------------------------------
x = x0;


% line search parameters
maxlsiter = params.lineSearchItnlim ; % Maximum number of linesearches that can happen
gradToll = params.gradToll ; % Gradient tolerance, I think - is used on a line below as a stopping criteria with respect to the 2-norm of dx
alpha = params.lineSearchAlpha;    beta = params.lineSearchBeta; % Line search criteria
t0 = params.lineSearchT0; % Starting value - default is 1
k = 0; % Counter
t = 0; % To be used in the first objective function down below

% compute g0  = grad(Phi(x))
g0 = wGradient(x,params);
% The wGradient is the overall change within the data x=x0
% This looks at the gradients in all senses of the gradient - it's looking
% at the object gradient [consistency of the data], the TV gradient and the
% XFM (Wavelet) gradient

dx = -g0; % The dx is the negative of the gradient calculated in the function
imshow(abs(dx)) % Show what the gradient actually looks like - Easy to see edges in this representation


% iterations
while (k > params.Itnlim) | (norm(dx(:)) < gradToll) 

% backtracking line-search

	% pre-calculate values, such that it would be cheap to compute the objective
	% many times for efficient line-search
	[FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx] = preobjective(x, dx, params);
        %The preobjective function is used to "precaclulate" the transforms
        %This calculates the transforms to be used after - I don't quite
        %understand the purpose of this, however.
    
    
	f0 = objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx,x,dx, t, params);
	t = t0;
        [f1, ERRobj, RMSerr]  =  objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx,x,dx, t, params);
	
	lsiter = 0; % Begin the counter

	while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:)))^2 & (lsiter<maxlsiter)
		lsiter = lsiter + 1;
		t = t * beta;
		[f1, ERRobj, RMSerr]  =  objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx,x,dx, t, params);
	end

	if lsiter == maxlsiter
		disp('Reached max line search,.... not so good... might have a bug in operators. exiting... ');
		return;
	end

	% control the number of line searches by adapting the initial step search
	if lsiter > 2
		t0 = t0 * beta;
	end 
	
	if lsiter<1
		t0 = t0 / beta;
	end

	x = (x + t*dx); % Apply the changes that we just calculated for! This is the part of the code that does the iterative CS

	%--------- uncomment for debug purposes ------------------------	
	%disp(sprintf('%d   , obj: %f, RMS: %f, L-S: %d', k,f1,RMSerr,lsiter));

	%---------------------------------------------------------------
	
    %conjugate gradient calculation
    
	g1 = wGradient(x,params); % Recalculate the gradient with the new information added onto x
	bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps); % bk is the ratio between the sum(g1(i)^2)/sum(g0(i)^2)
	g0 = g1; % The new gradient becomes the reference gradient
	dx =  - g1 + bk* dx; % The new dx is now the -new grad plus this ratio
	k = k + 1; % Counter increment
	
% 	%TODO: need to "think" of a "better" stopping criteria ;-)
% 	if (k > params.Itnlim) | (norm(dx(:)) < gradToll) 
% 		break;
% 	end

end


return;







