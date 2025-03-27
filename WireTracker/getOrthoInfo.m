%% get furthest point vector and final image of orthographic projection
% Jade Lariviere | last modified Mar. 24, 2025

function [direction_XYZ,stepSize,Image] = getOrthoInfo(sliceData,imageData,SearchType)
% function takes an image and, after performing image processing specific
% to this function's purpose using sliceData & imageData, returns the 
% direction [x,y,z] to the furthest point in the orthographic image in the
% best ROI and scaling stepSize. also returns the final image.

mult    = 3; % multiplier for binarization strictness; higher = stricter

img_base = imageData.base; img_mask = imageData.mask; 
size_img = imageData.dims; offset = imageData.offset;
axis_U = sliceData.axes.U; axis_V = sliceData.axes.V; % UV -> XYZ
% P is origin of UV-axis (0,0), which was offset in makeImage()
P_uv = [imageData.offset imageData.offset];

% finalize image ==========================================================
% apply mask over base (filled) image
img_final = img_mask.*(img_base~=0);
% isolate ROI which includes P(u,v)!
ROI_all = imbinarize(img_final,graythresh(img_final)*mult); % high thresh
labels = bwlabeln(ROI_all,8); % label ROIs
try
    label_ROI = labels(round(P_uv(1)),round(P_uv(2))); % ROI with P_uv?
catch
    warning('P(u,v) rounded outside of image bounds. Consider increasing ImageDims.');
    label_ROI = labels(min(P_uv(1),size_img(1)),min(P_uv(2),size_img(2)));
end

% find pixel with highest out of plane separation: makes your next vector!
ROI_best = img_final.*(labels==label_ROI);
[best_separation, best_uv, ~] = findMax(ROI_best,P_uv,7,50,SearchType);

% with plane equation + separation, point -> vector [x,y,z]
best_xyz = (best_uv(1)-offset)*axis_U + (best_uv(2)-offset)*axis_V;
best_dir = best_xyz + best_separation*sliceData.N; % direction vector!

% return furthest direction, step size & final image ======================
direction_XYZ = best_dir/norm(best_dir); % normalized direction
stepSize = abs(log(norm(best_dir))); % roughly scale step based on distance
Image.combined = img_final;
    Image.ROI = ROI_best;

% plots for debugging =====================================================
%{
labels_c = label2rgb(labels,'hsv','k','noshuffle'); % add color to plot
max_val = 2*max(img_final,[],'all');
valRange = [0 max_val*1.1]; % intensity scaling for plotting
ROI_best(round(best_uv(1)),round(best_uv(2))) = max_val; % highlight pixel

figure(11); % transpose images for plotting (MATLAB hates being normal)
    subplot(2,3,1); imshow(img_base',valRange); set(gca,"YDir","normal");
        xlabel('u'); ylabel('v'); title('Raw Ortho Projection');
    subplot(2,3,2); imshow(img_mask',valRange); set(gca,"YDir","normal");
        title('Image Mask (Pixel Fill)');
    subplot(2,3,3); imshow(img_final',valRange);
        set(gca,"YDir","normal"); title('Masked Base');
    subplot(2,3,4); imshow(ROI_all'); set(gca,"YDir","normal");
        title('Binarized')
    subplot(2,3,5); imshow(permute(labels_c,[2 1 3]));
        set(gca,"YDir","normal"); title('ROIs');
    subplot(2,3,6); imshow(ROI_best',valRange); set(gca,"YDir","normal");
        title('Best ROI','Furthest Distance Highlighted');
    colormap("pink");
%}
end
