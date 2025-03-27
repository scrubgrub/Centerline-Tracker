%% helper function to make 2D images from projected points
% Jade Lariviere | last modified Mar. 5, 2025

function [Image] = makeImage(sliceData,imgSize,ImageType)
% function takes sliceData with sorted coordinates [u,v], intensities, and
% other axis information to produce a 2D base image and image mask with 
% dimension imgSize. images are compatible with the Image Analysis Toolbox.

u = sliceData.pointsUV(:,1); v = sliceData.pointsUV(:,2); % pt coordinates
numPts = numel(u);
switch ImageType % select data represented by pixel intensity
    case 'Slice'; value = sliceData.intensity;
    case 'Ortho'; value = sliceData.separation;
    otherwise
        error('Input argument %s for ImageType not recognized.',ImageType);
end

% ensure ALL point coordinates are positive. coordinates = pixels
offset = abs(min(min(min(u),min(v)),0)) + 5; % apply positive offset
    u = u + offset; v = v + offset;
    
% create images ===========================================================
% logic prioritizes closer pixels than further ones due to sorted order
img_base = zeros(imgSize); % create blank base image
img_mask = img_base; % create identical img size mask over image

for i = 1:numPts
    % image base: rounds (u,v) to nearest integer as image pixel
    img_base(round(u(i)),round(v(i))) = value(i);

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
end