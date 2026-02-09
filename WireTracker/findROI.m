%% determine ROI closest to given point on image
% Jade Lariviere | last modified Aug. 12, 2025

function [LabelImage,numLabel,LabelBest] = findROI(Image,P_uv,connectivity)
% this helper function determines the closest ROI from a binarized Image
% given a particular pixel P [u,v] and pixel connectivity setting and 
% returns the labeled Image, number of labels, and the label of the ROI 
% closest to P.  ==========================================================

img_dims = size(Image);
[LabelImage,numLabel] = bwlabeln(Image,connectivity); % labeled ROIs of Image

% try 1: does coordinate fall in a ROI? -----------------------------------
LabelBest = LabelImage(min(ceil(P_uv(1)),img_dims(1)), ...
    (min(ceil(P_uv(2)),img_dims(2))));
if LabelBest ~= 0, return; end

% try 2: find closest ROI likely containing P_uv w/extrema ----------------
ROI_props = regionprops(LabelImage,'Extrema');
    ROI_props = cat(3,ROI_props.Extrema);
    size_ROIprops = size(ROI_props);
    if  numLabel > 1 % select closest ROI (apparent branching)
        % vector math: calculate shortest distance from P to ROI extrema
        UV = repmat(P_uv,[size_ROIprops(1) 1 size_ROIprops(3)]);
        [~,best_extrema] = min(vecnorm(ROI_props-UV,2,2),[],"all");
        [~,~,LabelBest] = ind2sub(size_ROIprops,best_extrema);
    elseif numLabel == 1 % presume only found ROI is correct
        LabelBest = 1; 
    else, LabelBest = 0; % no ROIs found on slice...
    end
% plots for debugging =====================================================
%{
labels_c = label2rgb(LabelImage,'hsv','k','noshuffle'); % add color to plot
final_img = LabelImage; final_img(final_img~=LabelBest) = 0;

figure(20); 
    subplot(1,3,1); imshow(Image'); set(gca,"YDir","normal");
        xlabel('u'); ylabel('v'); title('Input Image');
    subplot(1,3,2); imshow(permute(labels_c,[2 1 3])); set(gca,"YDir","normal");
        xlabel('u'); ylabel('v'); title('Identified ROIs');
    subplot(1,3,3); imshow(final_img'); set(gca,"YDir","normal");
        xlabel('u'); ylabel('v'); title('Best ROI Near P_uv');
%}
end