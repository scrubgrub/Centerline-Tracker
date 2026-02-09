%% helper function to make 2D images from projected points
% Jade Lariviere | last modified Jun. 23, 2025

function [Image] = makeImage(sliceData,imgSize,ImageType)
% function takes sliceData with sorted coordinates [u,v], intensities, and
% other axis information to produce a 2D base image and image mask with 
% dimension imgSize. images are compatible with the Image Analysis Toolbox.
% =========================================================================

u = sliceData.pointsUV(:,1); v = sliceData.pointsUV(:,2); % pt coordinates
switch ImageType % select data represented by pixel intensity
    case 'Slice'; value = sliceData.intensity;
    case 'Ortho'; value = sliceData.separation;
    otherwise
        error('Input argument %s for ImageType not recognized.',ImageType);
end

% prep data for image creation ============================================
% offset the projection so original origin is in the middle of image space
offset = [ceil(imgSize(1)/2) ceil(imgSize(2)/2)];
    u = u + offset(1); v = v + offset(2);
% discard points outside the image space boundaries
[pts_in] = inpolygon(u,v, ...
    [1 1 imgSize(1) imgSize(1)],[1 imgSize(2) imgSize(2) 1]);
    u(~pts_in) = []; v(~pts_in) = []; value(~pts_in) = []; % delete data
numPts = numel(u); % update number of elements tracked from now on

% create images ===========================================================
% getProjection logic prioritizes closer pixels due to sorted order
img_base = zeros(imgSize); % create blank base image
img_mask = img_base; % create identical img size to mask over image

for i = 1:numPts
    % image base: rounds (u,v) to nearest integer as image pixel
    img_base(min(ceil(u(i)),imgSize(1)), ...
        min(ceil(v(i)),imgSize(2))) = value(i);

    % image mask: ensures no background pixels infiltrate by spilling into
    % adjacent pixels (top, bottom, left and right)!
    img_mask(floor(u(i)),floor(v(i))) = value(i);
    img_mask(ceil(u(i)),floor(v(i))) = value(i);
    img_mask(floor(u(i)),ceil(v(i))) = value(i);
    img_mask(ceil(u(i)),ceil(v(i))) = value(i);
end
% fill any holes in image due to rounding
img_baseFilled = imfill(img_base,"holes");

% crop images to original dimensions, return Image structure ==============
Image.base = img_baseFilled(1:imgSize(1),1:imgSize(2));
    Image.mask = img_mask(1:imgSize(1),1:imgSize(2));
    Image.offset = offset; % store info of positive offset of coordinates
    Image.dims = size(Image.base); % store image dims

% plots for debugging =====================================================
%{
max_val = max(value,[],'all');
img_final = Image.mask.*(Image.base~=0);
valRange = [0 max_val*1.1]; % intensity scaling for plotting

figure(11); % transpose images for plotting (MATLAB hates being normal)
    subplot(1,2,1); imshow(img_final,valRange); set(gca,"YDir","normal");
        xlabel('u'); ylabel('v'); title('Final Raw Image');
    subplot(1,2,2); histogram(img_final,50, "Normalization","percentage");
        xlabel('Value Distribution'); ylabel('Proportion');
%}
end
