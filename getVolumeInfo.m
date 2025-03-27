%% get local max "centroid" coordinates and final image of cross-section
% Jade Lariviere | last modified Mar. 24, 2025

function [max_XYZ,Image] = getVolumeInfo(sliceData,imageData,SearchType)
% function takes an image and, after performing image processing specific
% to this function's purpose using sliceData & imageData, returns the 
% coordinates [x,y,z] of the local maximum intensity closest to the
% starting point (point P in UV space). also returns the final image.

mult    = 0.9; % multiplier for binarization strictness; higher = stricter

img_base = imageData.base; img_mask = imageData.mask;
img_dims = imageData.dims;
axis_U = sliceData.axes.U; axis_V = sliceData.axes.V; % UV -> XYZ
% P is origin of UV-axis (0,0), which was offset in makeImage()
offset = imageData.offset; P_uv = [offset offset];

% finalize image ==========================================================
% apply mask over base image, prioritizing larger intensities
img_final = img_base.*(img_base >= img_mask) + img_mask.*(img_mask > img_base);
    img_final = img_final.*(img_base~=0); % mask constrained to base

% isolate ROI which includes P_uv =========================================
ROI_all = imbinarize(img_final,graythresh(img_final)*mult);
labels = bwlabeln(ROI_all,8); % label ROIs
try
    label_ROI = labels(round(P_uv(1)),round(P_uv(2))); % ROI with P_uv?
catch
    warning('P(u,v) rounded outside of image bounds. Consider increasing ImageDims.');
    label_ROI = labels(min(P_uv(1),img_dims(1)),min(P_uv(2),img_dims(2)));
end
ROI_best = img_final.*(labels==label_ROI);

% find local max intensity ================================================
[~,max_uv,cMap] = findMax(ROI_best,P_uv,5,20,SearchType);

% convert UV -> XYZ space, return XYZ & final Image =======================
% remove offset
max_XYZ = (max_uv(1)-offset)*axis_U + (max_uv(2)-offset)*axis_V + sliceData.P;
Image.combined = img_final;
    Image.ROI = ROI_best;

% plots for debugging =====================================================
%{
labels_c = label2rgb(labels,'hsv','k','noshuffle'); % add color to plot
max_val = 2*max(img_final,[],'all');
valRange = [0 max_val*1.1]; % intensity scaling for plotting
ROI_best(round(max_uv(1)),round(max_uv(2))) = max_val; % highlight found pixel

figure(12); % transpose images for plotting (MATLAB hates being normal)
    subplot(2,3,1); imshow(img_base',valRange); set(gca,"YDir","normal");
        xlabel('u'); ylabel('v'); title('Raw Cross-section');
    subplot(2,3,2); imshow(img_mask',valRange); set(gca,"YDir","normal");
        title('Image Mask (Pixel Fill)');
    subplot(2,3,3); imshow(img_final',valRange);
        set(gca,"YDir","normal"); title('Masked Base (Max)');
    subplot(2,3,4); imshow(permute(labels_c,[2 1 3])); 
        set(gca,"YDir","normal"); title('ROIs');
    subplot(2,3,5); imshow(permute(cMap,[2 1 3])); 
        set(gca,"YDir","normal"); title('Derivative Colormap');
    subplot(2,3,6); imshow(ROI_best',valRange); set(gca,"YDir","normal");
        title('Best ROI','Centroid Highlighted');
%}
end