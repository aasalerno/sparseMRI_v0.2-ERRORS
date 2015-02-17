function [minIntrVec,stat,actpctg] = genSampling(pdf,iter,tol)

%[mask,stat,N] = genSampling(pdf,iter,tol)
%
% a monte-carlo algorithm to generate a sampling pattern with 
% minimum peak interference. The number of samples will be
% sum(pdf) +- tol
%
%	pdf - probability density function to choose samples from
%	iter - number of tries
%	tol  - the deviation from the desired number of samples in samples
%
% returns:
%	mask - sampling pattern
%	stat - vector of min interferences measured each try
%	actpctg    - actual undersampling factor
%
%	(c) Michael Lustig 2007

h = waitbar(0);

pdf(find(pdf>1)) = 1; % Double check to make sure that there's nothing greater than 1
K = sum(pdf(:)); % Equvalent to N from genPDF - number of elements = sum of PDF

minIntr = 1e99; % Absurdly large value, as a preemptive
minIntrVec = zeros(size(pdf)); % Preallocate a vector for the data

for n=1:iter
	tmp = zeros(size(pdf)); % preallocate
	while abs(sum(tmp(:)) - K) > tol
		tmp = rand(size(pdf))<pdf; % creates a map of values
        % If it is LOWER than the value created in the PDF map, then it's
        % like that one has been chosen. Our K value is the amount that
        % we're allowed to have, so we are trying to build a map of a small
        % enough amount of allowed values, within tolerance (i.e. not too
        % few and not too many)
	end
	
	TMP = ifft2(tmp./pdf);
	if max(abs(TMP(2:end))) < minIntr % First value is the DC offset - either way this is always true
		minIntr = max(abs(TMP(2:end)));
		minIntrVec = tmp;
	end
	stat(n) = max(abs(TMP(2:end))); % This tells us the max value for each iteration
	waitbar(n/iter,h);
end

actpctg = sum(minIntrVec(:))/prod(size(minIntrVec)); % average value

close(h); % closes the waitbar


