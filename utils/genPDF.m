function [pdf,val] = genPDF(imSize,p,pctg,distType,radius,disp)

%[pdf,val] = genPDF(imSize,p,pctg [,distType,radius,disp])
%
%	generates a pdf for a 1d or 2d random sampling pattern
%	with polynomial variable density sampling
%
%	Input:
%		imSize - size of matrix or vector - this is the dataset size
%		p - power of polynomial - how is the falloff occurring
%		pctg - partial sampling factor e.g. 0.5 for half
%		distType - 1 or 2 for L1 or L2 distance measure
%		radius - radius of fully sampled center - THIS IS A DECIMAL OF THE
%		WHOLE - i.e cannot be greater than 1
%		disp - display output
%
%	Output:
%		pdf - the pdf
%		val - min sampling density
%
% 
%	Example:
%	[pdf,val] = genPDF([128,128],2,0.5,2,0,1);
%
%	(c) Michael Lustig 2007

if nargin < 4
	distType = 2; %L2 distance measuring
end

if nargin < 5
	radius = 0; % Radius of full sampling
end

if nargin < 6
	disp = 0; % Display output
end


minval=0;
maxval=1;
val = 0.5;

if length(imSize)==1
	imSize = [imSize,1]; %Make a 1D image 2D if required
end

sx = imSize(1); 
sy = imSize(2);
PCTG = floor(pctg*sx*sy); %Dictates the number of points that need to be chosen


if sum(imSize==1)==0  % 2D - i.e. if something is imSize==1, it gives a 1
	[x,y] = meshgrid(linspace(-1,1,sy),linspace(-1,1,sx)); % Ensures that the zero is in the middle
	switch distType % Looking at which L1 or L2 is wanted
		case 1 % Does this if the L1 system is wanted
			r = max(abs(x),abs(y));
        otherwise % Does this if the L2 is wanted
			r = sqrt(x.^2+y.^2); 
			r = r/max(abs(r(:))); % This makes a small centre and large outside
	end

else %1d
	r = abs(linspace(-1,1,max(sx,sy)));
end




idx = find(r<radius); % Looks for values that satisfy the condition that
                      % their value relative to the centre is less than the
                      % value of the radius. It's clever in the fact that
                      % it looks at the values as a distribution with
                      % respect to the centre, as a distance function

pdf = (1-r).^p; % Here, the higher the power, the better we are, as it will decrease the values lower
pdf(idx) = 1;
if floor(sum(pdf(:))) > PCTG
	error('infeasible without undersampling dc, increase p');
end

% begin bisection
while(1)
	val = minval/2 + maxval/2; 
	pdf = (1-r).^p + val; % Seems to arbitrarily add a Prob Density
    pdf(find(pdf>1)) = 1; % Reset to viable values - this arbitrary add on may also bring other values high enough
    pdf(idx)=1;
	N = floor(sum(pdf(:)));
    % Here we have our sanity checking... We want the sum of the pdf to be
    % exactly the number of points we need, making us know a bit about how
    % the data will be sampled, I think...
	if N > PCTG % infeasible
		maxval=val;
	end
	if N < PCTG % feasible, but not optimal
		minval=val;
	end
	if N==PCTG % optimal - ONE optimization I'd do here is some tolerance bounds
		break;
	end
end

if disp
	figure, 
	subplot(211), imshow(pdf)
	if sum(imSize==1)==0
		subplot(212), plot(pdf(end/2+1,:));
	else
		subplot(212), plot(pdf);
	end
end






