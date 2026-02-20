%% get furthest point vector and final image of orthographic projection
% Jade Lariviere | last modified Jan. 21, 2026

function [DirectionXYZ,stepSize,Image] = getOrthoInfo(sliceData,imageData,SearchType)
% function performs image processing using sliceData & imageData to make
% a "camera" image representation and determine travel direction 
% based on image pixel depth. returns the XYZ direction vector to the 
% selected point in the orthographic camera, weighted variable stepSize
% based on how far said point is from current position P, and the Image
% structure containing final images and additional information. argument 
% SearchType sets the findMax function method. ============================

mult    = 2.6;  % multiplier for binarization strictness; higher = stricter
connect = 4;    % specify 2D connectivity of bwlabeln (strict)
maxFlag = 7;    % for findMax; max number of times threshold flag triggers
maxIt   = 50;   % for findMax; max number of search iterations

img_base = imageData.base; img_mask = imageData.mask; 
img_dims = imageData.dims; 
% P is origin of UV-axis (0,0), which was offset in makeImage()
P_uv = imageData.offset; 

% finalize image ==========================================================
% apply mask over base (filled) image
img_final = img_mask.*(img_base~=0);
min(img_final(img_final <= 0),[],"all");
% isolate ROI which includes P(u,v)!
ROI_all = imbinarize(img_final,graythresh(img_final)*mult); % high thresh
[labels,~,label_ROI] = findROI(ROI_all,P_uv,connect);

% find pixel with highest out of plane separation: makes your next vector!
ROI_best = img_final.*(labels==label_ROI);
[best_separation,best_uv,~] = findMax(ROI_best,P_uv,maxFlag,maxIt,SearchType);

% with plane equation + separation, point -> vector [x,y,z]
best_dir = UVtoXYZ(best_uv,sliceData.axes,best_separation*sliceData.N,P_uv);
    if best_dir == [0,0,0], best_dir = sliceData.N; end

% return furthest direction, step size & final image ======================
DirectionXYZ = best_dir/norm(best_dir); % normalized direction
stepSize = max(abs(log(norm(best_dir))),1); % scale step with magnitude
Image.all = img_final; % output Image info
    Image.ROI = ROI_best;
    Image.range = [0 max(img_final,[],"all")];

% plots for debugging =====================================================
%{
labels_c = label2rgb(labels,'hsv','k','noshuffle'); % add color to plot
max_val = 2*max(img_final,[],'all');
valRange = [0 max_val*1.1]; % intensity scaling for plotting
ROI_best(min(ceil(best_uv(1)),img_dims(1)), ...
    min(ceil(best_uv(2)),img_dims(2))) = max_val; % highlight pixel

figure(12); % transpose images for plotting (MATLAB hates being normal)
    subplot(2,3,1); imshow(img_base',valRange); set(gca,"YDir","normal");
        xlabel('u'); ylabel('v'); title('Raw Image');
    subplot(2,3,2); imshow(img_mask',valRange); set(gca,"YDir","normal");
        title('Neighbor Pixels Flooded');
    subplot(2,3,3); imshow(img_final',valRange);
        set(gca,"YDir","normal"); title('Flood Masked to Raw');
    subplot(2,3,4); imshow(ROI_all'); set(gca,"YDir","normal");
        title('Binarized Image');
    subplot(2,3,5); imshow(permute(labels_c,[2 1 3]));
        set(gca,"YDir","normal"); title('Found ROIs');
    subplot(2,3,6); imshow(ROI_best',valRange); set(gca,"YDir","normal");
        title({'Best Direction','Furthest Pixel Highlighted'});
    colormap("pink");
%}
end
