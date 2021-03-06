function [img,fil] = testMap(file,sty,sl,loc,dir)
%
% testMap.m                         02/18/15
%
% This function is built so that we can do a recon on the data in a
% specific orientation of our choosing and go from there.
% It assumes that the RO direction is 3.
%
% Make sure that we have all the required information for the strcmp

if nargin < 4; loc = 1; end
if nargin < 5; dir = []; end

if strcmp(sty,'par') || strcmp(sty,'per')
    if isempty(dir)
        error('You must specify a direction if using per or par');
    end
    if strcmp(sty,'per')
        dir = [-dir(2) dir(1)];
    end
end

if isa(file,'cell')
    data = mincmap(file{1})+1i*mincmap(file{2});
    data = fftshift(fft(data,[],loc),loc);
    % Process it for powers of two
    dims = size(data);
    
    dims2 = 2.^(nextpow2(dims)-1);
    
    for i = 1:3
        if dims(i) == 2*dims2(i)
            dims2(i) = dims(i);
        end
        if mod(dims(i),2) ~= 0
            dims(i) = dims(i)+1;
            if i == 1
                data = data(1:end-1,:,:);
            elseif i == 2
                data = data(:,1:end-1,:);
            elseif i == 3
                data = data(:,:,1:end-1);
            end
        end
    end
    
    res = ceil((dims - dims2)/2); % The residual that we need to see on either side in order to compress the data to 2^n
    
    dataRes = data(res(1)+1:end-res(1),res(2)+1:end-res(2),res(3)+1:end-res(3));
    if loc == 1
        dataUse = squeeze(dataRes(round(sl/dims(1)*dims2(1)),:,:));
    elseif loc == 2
        dataUse = squeeze(dataRes(:,round(sl/dims(2)*dims2(2)),:));
    elseif loc == 3
        dataUse = squeeze(dataRes(:,:,round(sl/dims(3)*dims2(3))));
    end
    if isa(dataUse,'uint16')
        dataUse = double(dataUse);
    end
    
    img = fftshift(fft2(dataUse));
    
%     % Build the basis of the filter
    if strcmp(sty,'par') || strcmp(sty,'per')
        slp = dir(2)/dir(1);
        fil = linefilt(dataUse,slp,sampFac/2);
    else
        fil = zeros(size(dataUse));
    end
%     elseif strcmp(sty,'circ')
%         fil = circfilt(dataUse,sampFac/2);
%     elseif strcmp(sty,'lores') || strcmp(sty,'sq')
%         fil = sqfilt(dataUse,sampFac/2);
%     end
%     
%     
%     sz = size(fil);
%     % Add the 1/r component
%     alim = 1/sqrt(2);
%     [x,y] = meshgrid(linspace(-alim,alim,sz(1)),linspace(-alim,alim,sz(2)));
%     r = sqrt(x.^2+y.^2);
%     fil2 = (ones(sz)-fil).*rand(sz(1),sz(2)).*(1-r/max(r(:)));
%     filTot = fil + fil2; % Build the total filter PDF
%     pdf = (1-r).*(ones(sz)-fil)+fil;
%     for i = 1:sz(1)
%         for j = 1:sz(2)
%             if filTot(i,j) < (0.75*sampFac)
%                 filTot(i,j) = 1;
%             else
%                 filTot(i,j) = 0;
%             end
%         end
%     end
%     % We know that if two files are given, one is the real, the other
%     % is the imaginary
%     img = dataUse.*filTot;
    
    
else
%     % Read in the data
%     data = mincread(file,'image');
%     
%     % Process it for powers of two
%     dims = size(data);
%     
%     dims2 = 2.^(nextpow2(dims)-1);
%     
%     for i = 1:3
%         if dims(i) == 2*dims2(i)
%             dims2(i) = dims(i);
%         end
%         if mod(dims(i),2) ~= 0
%             dims(i) = dims(i)+1;
%             if i == 1
%                 data = data(1:end-1,:,:);
%             elseif i == 2
%                 data = data(:,1:end-1,:);
%             elseif i == 3
%                 data = data(:,1:end-1,:);
%             end
%         end
%     end
%     
%     res = ceil((dims - dims2)/2); % The residual that we need to see on either side in order to compress the data to 2^n
%     
%     dataRes = data(res(1)+1:end-res(1),res(2)+1:end-res(2),res(3)+1:end-res(3));
%     if loc == 1
%         dataUse = squeeze(dataRes(round(sl/dims(1)*dims2(1)),:,:));
%     elseif loc == 2
%         dataUse = squeeze(dataRes(:,round(sl/dims(2)*dims2(2)),:));
%     elseif loc == 3
%         dataUse = squeeze(dataRes(:,:,round(sl/dims(3)*dims2(3))));
%     end
%     if isa(dataUse,'uint16')
%         dataUse = double(dataUse);
%     end
%     
%     
%     % Build the basis of the filter
%     if strcmp(sty,'par') || strcmp(sty,'per')
%         slp = dir(2)/dir(1);
%         fil = linefilt(dataUse,slp,sampFac/2);
%     elseif strcmp(sty,'circ')
%         fil = circfilt(dataUse,sampFac/2);
%     elseif strcmp(sty,'lores') || strcmp(sty,'sq')
%         fil = sqfilt(dataUse,sampFac/2);
%     end
%     
%     sz = size(fil);
%     % Add the 1/r component
%     alim = 1/sqrt(2);
%     [x,y] = meshgrid(linspace(-alim,alim,sz(1)),linspace(-alim,alim,sz(2)));
%     r = sqrt(x.^2+y.^2);
%     fil2 = (ones(sz)-fil).*rand(sz(1),sz(2)).*(1-r/max(r(:)));
%     filTot = fil + fil2; % Build the total filter PDF
%     for i = 1:sz(1)
%         for j = 1:sz(2)
%             if filTot(i,j) < (0.75*sampFac)
%                 filTot(i,j) = 0;
%             else
%                 filTot(i,j) = 1;
%             end
%         end
%     end
%     
%     if ksp == 0
%         % Normalize to be safe
%         imshow(dataUse,[]);
%         dataUse = dataUse./max(dataUse(:));
%         dataUse = fftshift(fft2(dataUse));
%     end
%     img = dataUse.*filTot;
%     
%     if 1 == 1
%         figure; imshow(abs(dataUse),[]);
%         figure; imshow(filTot,[]);
%     end
end