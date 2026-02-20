%% helper function to detect closest local maximum of 2D image
% Jade Lariviere | last modified Jan. 21, 2026

function [maxVal,maxUV,cMap] = findMax(Image,startPt,maxFlags,maxIt,SearchType)
% this function calculates the gradient of the input Image and iterates
% based on flagSensitivity or until MaxIt to find the closest local max. it
% returns the maximum value found and the associated coordinates (u,v). the
% SearchType determines what method to use to determine the maximum value
% in an image ROI, which can be customized. ===============================

scaleFactor = 5; % how much bigger to make upscaled img; 1 = same
interpMethod = 'bicubic'; % interp method for imresize

% make higher pixel image via interp, scale startPt
[newImage,newImg_dim] = scaleImg(Image,scaleFactor,interpMethod);
newStartPt = startPt.*scaleFactor;
img_dim = size(Image); % store original image dims

% calculate image gradient from best ROI ----------------------------------
[fv,fu] = gradient(newImage); % swap fv/fu due to gradient() definition

% RGB color map for viewing image derivatives
cMap(1:newImg_dim(1),1:newImg_dim(2),1:3) = 0;
    c_fu = fu; c_fv = fv;
    c_fu(fu<0) = fu(fu<0)./abs(min(fu.*2,[],'all'));
    c_fu(fu>0) = fu(fu>0)./max(fu.*2,[],'all');
    c_fv(fv<0) = fv(fv<0)./abs(min(fv.*2,[],'all'));
    c_fv(fv>0) = fv(fv>0)./max(fv.*2,[],'all');
    cMap(:,:,2) = c_fu + 0.5; cMap(:,:,3) = c_fv + 0.5;

% Search Type =============================================================
switch SearchType % determine which max determination method to use
case 'Absolute'
% grab max value ----------------------------------------------------------
[maxVal,idx_max] = max(newImage,[],'all');
[max_u, max_v] = ind2sub(newImg_dim,idx_max);
    maxUV = [max_u, max_v]./scaleFactor; % scale back to Image dims!
case 'Local'
% search loop -------------------------------------------------------------
c_new = newStartPt; % starting coordinates (P = origin)
try value_new = newImage(ceil(c_new(1)),ceil(c_new(2))); % pixel value
catch % clamp coordinates to image dimensions
    warning('P(u,v) indexed past image dimensions. Consider increasing ImageDims.');
    c_new = [min(c_new(1),newImg_dim(1)), min(c_new(2),newImg_dim(2))]; % clamp
    value_new = newImage(c_new(1),c_new(2));
end

flagCount = 0; flagThresh = maxFlags; i = 0;
while (flagCount < flagThresh) && (i < maxIt)
    % update pixel (centroid) & intensity value from previous loop
    c = c_new; value = value_new;
    % store direction of ascent
    ascent = [fu(min(ceil(c(1)),newImg_dim(1)),min(ceil(c(2)),newImg_dim(2))), ...
        fv(min(ceil(c(1)),newImg_dim(1)),min(ceil(c(2)),newImg_dim(2)))];
    ascent_dir = ascent/vecnorm(ascent,2,2); % normalize

    if isnan(ascent_dir) % what if ascent of zero at c?
        % check neighboring pixels for valid ascent, average
        neighbors_ascent = zeros(3,3,2); % 3x3 grid for 1) fu & 2) fv
        neighbors_ascent(:,:,1) = fu(min(ceil(c(1)-1:c(1)+1),newImg_dim(1)), ...
            min(ceil(c(2)-1:c(2)+1),newImg_dim(2)));
        neighbors_ascent(:,:,2) = fv(min(ceil(c(1)-1:c(1)+1),newImg_dim(2)), ...
            min(ceil(c(2)-1:c(2)+1),newImg_dim(2)));
        neighbors_ascent(isnan(neighbors_ascent)) = 0; % NaN failsafe

        ascent = squeeze(mean(neighbors_ascent,[1 2]))'; % average ascent
            ascent_dir = ascent_dir/vecnorm(ascent,2,2); % normalize
        if isnan(ascent_dir), ascent_dir = [0 0]; end % final failsafe
    end
    % step to new pixel (log scaling)
    c_new = c + ascent_dir.*(abs(log(vecnorm(ascent,2,2)+5))./4);
        % clamp c to within image bounds to prevent index errors!
        c_new = [max(min(c_new(1),newImg_dim(1)),1) max(min(c_new(2),newImg_dim(2)),1)];
    % what is the value at the pixel where c is?
    value_new = newImage(min(ceil(c_new(1)),newImg_dim(1)), ...
        min(ceil(c_new(2)),newImg_dim(2)));
    % compare new_value to value; always seek higher value
    if value_new <= value
        flagCount = flagCount + 1; % settled or overshot, track it
    end
    i = i + 1; % increment counter; repeat while loop conditions are true
end
% return final max coordinates + value ====================================
% scale back to original Image dimensions + use intensity in original
maxUV = c_new./scaleFactor;
maxVal = Image(min(ceil(maxUV(1)),img_dim(1)), ...
    min(ceil(maxUV(2)),img_dim(2))); % 
otherwise
error('Input argument %s for SearchType not recognized.',SearchType);
end

% plotting for debugging ==================================================
%{
img = Image;
max_val = 2*max(img,[],'all');
valRange = [0 max_val*1.1]; % intensity scaling for plotting
derivRange = [min([fu fv],[],"all")*1.1, max([fu fv],[],"all")*1.1];
% highlight found pixel
img(min(ceil(maxUV(1)),img_dim(1)), ...
    min(ceil(maxUV(2)),img_dim(2))) = max_val;
    
figure(14); % transpose images for plotting (MATLAB hates being normal)
    subplot(2,2,1); imshow(img',valRange); set(gca,"YDir","normal");
        title('Original Image (POI Highlighted)');
    subplot(2,2,2); imshow(permute(cMap,[2 1 3])); set(gca,"YDir","normal");
        title('Derivative Colormap of Scaled Image');
    subplot(2,2,3); imshow(fu',derivRange); set(gca,"YDir","normal");
        xlabel('u'); ylabel('v'); title('F_u');
    subplot(2,2,4); imshow(fv',derivRange); set(gca,"YDir","normal");
        xlabel('u'); ylabel('v'); title('F_v');
%}
end
