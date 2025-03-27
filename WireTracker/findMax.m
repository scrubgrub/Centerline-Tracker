%% helper function to detect closest local maximum of 2D image
% Jade Lariviere | last modified Mar. 24, 2025

function [maxVal,maxUV,cMap] = findMax(Image,startPt,maxFlags,maxIt,SearchType)
% this function calculates the gradient of an input 2D image and iterates
% based on flagSensitivity or until MaxIt to find the closest local max. it
% returns the maximum value found and the associated coordinates (u,v). the
% SearchType determines what method to use to determine the maximum value
% in an image ROI, which is either absolute or local search.

img_dims = size(Image); % grab image dimensions info

% calculate image gradient from best ROI ----------------------------------
[fv,fu] = gradient(Image); % swap fv/fu due to gradient() definition

% RGB color map for viewing image derivatives
cMap(1:img_dims(1),1:img_dims(2),1:3) = 0;
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
[maxVal,idx_max] = max(Image,[],'all');
[max_u, max_v] = ind2sub(img_dims,idx_max);
    maxUV = [max_u, max_v];
case 'Local'
% search loop -------------------------------------------------------------
c_new = startPt; % starting coordinates (P = origin)
try value_new = Image(round(c_new(1)),round(c_new(2))); % pixel value
catch % clamp coordinates to image dimensions
    warning('P(u,v) clamped to image dimensions. Consider increasing ImageDims.');
    c_new = [min(c_new(1),img_dims(1)), min(c_new(2),img_dims(2))]; % clamp
    value_new = Image(c_new(1),c_new(2));
end

flagCount = 0; flagThresh = maxFlags; i = 0;
while (flagCount < flagThresh) && (i < maxIt)
    % update pixel (centroid) & intensity value from previous loop
    c = c_new;
    value = value_new;
    
    ascent = [fu(round(c(1)),round(c(2))), ...
        fv(round(c(1)),round(c(2)))]; % store closest direction of ascent
    ascent_dir = ascent/vecnorm(ascent,2,2);
        if isnan(ascent_dir); ascent_dir = [0 0]; end % prevent NaNs
    % step to new pixel/centroid (logs cale clamped), get intensity value
    c_new = c + ascent_dir.*(abs(log(vecnorm(ascent,2,2)+5))./4);
    value_new = Image(round(c_new(1)),round(c_new(2)));
    % compare new_value to value; always seek higher value
    if value_new <= value
        flagCount = flagCount + 1; % overshot max, internally track it
    end
    i = i + 1; % increment counter; repeat while loop conditions are true
end
% return final max value + coordinates ====================================
maxVal = value_new; maxUV = c_new;
otherwise
error('Input argument %s for SearchType not recognized.',SearchType);
end
% plotting for debugging ==================================================
%{
max_val = 2*max(Image,[],'all');
valRange = [0 max_val*1.1]; % intensity scaling for plotting
Image(round(maxUV(1)),round(maxUV(2))) = max_val; % highlight found pixel
    
figure(13); % transpose images for plotting (MATLAB hates being normal)
    subplot(2,2,1); imshow(Image',valRange); set(gca,"YDir","normal");
        title('Image');
    subplot(2,2,2); imshow(permute(cMap,[2 1 3])); set(gca,"YDir","normal");
        title('Colormap of Image Derivative');
    subplot(2,2,3); imshow(fu',[-0.5 0.5]); set(gca,"YDir","normal");
        xlabel('u'); ylabel('v'); title('F_u');
    subplot(2,2,4); imshow(fv',[-0.5 0.5]); set(gca,"YDir","normal");
        xlabel('u'); ylabel('v'); title('F_v');
%}
end
