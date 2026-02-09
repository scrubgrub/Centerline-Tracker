%% get local max "centroid" coordinates and final image of cross-section
% Jade Lariviere | last modified Jan. 21, 2025

function [CenterXYZ,Image] = getVolumeInfo(sliceData,imageData,SearchType)
% function performs image processing using sliceData & imageData to make
% a cross section slice image and find a local maximum "centroid" based on 
% the image pixel intensity. returns the XYZ coordinates of the selected 
% centroid point in the cross section and the Image structure containing 
% final images and additional information. argument SearchType sets the
% findMax function method. ================================================

mult    = 0.85; % multiplier for binarization strictness; higher = stricter
radius  = 1;    % erosion structural element (disk) radius; larger = more
blur    = 2;    % standard deviation for gaussian filtering
connect = 8;    % specify 2D connectivity of bwlabeln
maxFlag = 5;    % for findMax; max number of times threshold flag triggers
maxIt   = 20;   % for findMax; max number of search iterations

img_base = imageData.base; img_mask = imageData.mask;
img_dims = imageData.dims;
% P is origin of UV-axis (0,0), which was offset in makeImage()
P_uv = imageData.offset; % P_uv essentially = offset

% finalize image ==========================================================
% apply mask over base image, prioritizing larger intensities
img_final = img_base.*(img_base >= img_mask) + img_mask.*(img_mask > img_base);
    img_final = img_final.*(img_base~=0); % mask constrained to base

% image processing ========================================================
img_final = double(bwmorph(img_final,'close'));
% erode & blur the cross section slice to create weighted volume
img_final = imerode(img_final,strel("disk",radius));
img_final = imgaussfilt(img_final, blur);

% isolate ROI which includes P_uv =========================================
ROI_all = imbinarize(img_final,graythresh(img_final)*mult);
[labels,~,label_ROI] = findROI(ROI_all,P_uv,connect); 
ROI_best = img_final.*(labels==label_ROI);

% find local max (centroid) ===============================================
[~,center_uv,cMap] = findMax(ROI_best,P_uv,maxFlag,maxIt,SearchType);

% convert UV -> XYZ space, return XYZ & final Image =======================
% remove offset (effectively P in UV coordinates)
CenterXYZ = UVtoXYZ(center_uv,sliceData.axes,sliceData.P,P_uv);
Image.all = img_final; % export Image information
    Image.ROI = ROI_best;
    Image.range = [0 max(img_final,[],"all")];

% plots for debugging =====================================================
%{
labels_c = label2rgb(labels,'hsv','k','noshuffle'); % add color to plot
max_val = max(img_final,[],'all');
valRange = [0 max_val*1.1]; % intensity scaling for plotting
ROI_best(min(ceil(center_uv(1)),img_dims(1)), ...
    min(ceil(center_uv(2)),img_dims(2))) = max_val*1.1; % highlight found pixel

figure(13); % transpose images for plotting (MATLAB hates being normal)
    subplot(2,3,1); imshow(img_base',valRange); set(gca,"YDir","normal");
        xlabel('u'); ylabel('v'); title('Raw Image');
    subplot(2,3,2); imshow(img_mask',valRange); set(gca,"YDir","normal");
        title('Neighbor Pixels Flooded');
    subplot(2,3,3); imshow(img_final',valRange);
        set(gca,"YDir","normal"); title({'Flood Masked to Raw', ...
            '+ Post-Processing'});
    subplot(2,3,4); imshow(permute(labels_c,[2 1 3])); 
        set(gca,"YDir","normal"); title('Found ROIs');
    subplot(2,3,5); imshow(permute(cMap,[2 1 3])); 
        set(gca,"YDir","normal"); title('Scaled Gradient Colormap');
    subplot(2,3,6); imshow(ROI_best',valRange); set(gca,"YDir","normal");
        title({'Best Centroid','Pointer Highlighted'});
%}
end
