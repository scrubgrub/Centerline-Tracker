%% rudimentary branch detection algorithm
% Jade Lariviere | last modified Jan. 22, 2026

function [branchVol, branchDir] = findBranch(branchVol,Slice,branchDir,I)
% function identifies what may be branching while tracking a path using
% image processing techniques on the cross section slice. the function
% takes the 2D Slice struct content, processes it to identify branches, and 
% updates the flagged pixels and the average direction vector pointing to 
% the pixels in XYZ space in the 3D matrices branchVol and branchDir, 
% respectively. ===========================================================

spacing     = 3; % grid density covering image ROI; larger = less points
flagMax     = 15; % max number of flag triggers before ending loop
maxIt       = 60; % max number of iterations for findMax
connect     = 8; % connectivity parameter

img = Slice.FinImage.ROI; % upscaled image from input argument
img_dims = Slice.Image.dims; % image dimensions
vol_dims = size(branchVol); % storage volume dimensions
P_uv = Slice.Image.offset; % to correct image offset
P = Slice.Data.P;

% prepare spaced grid of start points -> convert to array =================
grid_img = zeros(img_dims);
    grid_img((1:spacing:end),(1:spacing:end)) = 1;
    grid_img(~img) = 0; % mask grid to only points in image 
[grid_u, grid_v] = find(grid_img);

grid_pts = [grid_u grid_v]; % list of all points on grid

% run findMax to aggregate results into possible ROIs on Image ============
search_pts = nan(size(grid_pts));
search_vals = nan(length(search_pts),1);
search_img = zeros(img_dims);
for i = 1:size(grid_pts,1)
    % run findMax, record final coordinate + value
    [search_vals(i),search_pts(i,:)] = findMax(img,grid_pts(i,:), ...
        flagMax,maxIt,'Local');
    UV = [min(ceil(search_pts(i,1)),img_dims(1)), ...
        min(ceil(search_pts(i,2)),img_dims(2))]; % clamp to image indices

    % reinforcement logic: if pixel already = 1, flood adjacent pixels
    if search_img(UV(1),UV(2)) == 1
        % set adjacent pixels (3x3 area) true-ish
        try search_img(max(UV(1)-1,0):min(UV(1)+1,img_dims(1)), ...
                max(UV(2)-1,0):min(UV(2)+2,img_dims(2))) = 0.6;
            search_img(UV(1),UV(2)) = 1;
        catch
            warning('findBranches reinforcement logic index failure (i=%d).',I)
        end
    else % proceed like normal
        search_img(UV(1),UV(2)) = 1;
    end
end
% process final image -----------------------------------------------------
search_finalImg = search_img;
    search_finalImg = bwmorph(search_finalImg,'close');
    search_finalImg = bwmorph(search_finalImg,'fill');
    search_finalImg = bwmorph(search_finalImg,'open');

% find and remove closest ROI (where tracker should be) ===================
[labels,label_num,label_ROI] = findROI(search_finalImg,P_uv,connect);
    if label_num > 1 % remove closest ROI; exists apparent branching
        found_branches = search_finalImg.*(labels~=label_ROI);
    else, return; % if ROI = 1, that is the tracked wire; no branching
    end

% extract branching ROIs, convert UV -> XYZ, return updated volume ========
[vol_UV(:,1), vol_UV(:,2)] = find(found_branches);

% remove offset (effectively P in UV coordinates)
vol_XYZ = UVtoXYZ(vol_UV,Slice.Data.axes,P,P_uv);
    vol_XYZ = min(ceil(vol_XYZ),vol_dims,"includemissing");
separation_dir = vol_XYZ - P;
    separation_dir = separation_dir./vecnorm(separation_dir,2,2); % norm

for i = 1:size(vol_XYZ,1)
    % for each XYZ point of interest, set voxel to exist as memory
    xyz = vol_XYZ(i,:);     branchVol(xyz(1),xyz(2),xyz(3)) = 1;
    % try assigning neighboring voxels 1 to expand ROI area
    try branchVol(max(xyz(1)-1,1):min(xyz(1)+1,vol_dims(1)), ...
            max(xyz(2)-1,1):min(xyz(2)+1,vol_dims(2)), ...
            max(xyz(3)-1,1):min(xyz(3)+1,vol_dims(3))) = 1;
    catch
        warning('findBranches memory assignment of adjacent branching voxels failed (i=%d).',I);
    end
    % for each XYZ POI, store direction from camera view
    branchDir{vol_XYZ(i,1),vol_XYZ(i,2),vol_XYZ(i,3)} = separation_dir(i,:);
end
% plots for debugging =====================================================
%{
labels_c = label2rgb(labels,'hsv','k','noshuffle'); % add color
valRange = [0 max(img,[],'all')*1.1]; % intensity scaling for plotting

figure(15); % transpose images for plotting (MATLAB hates being normal)
    subplot(1,6,1); imshow(img',valRange); set(gca,"YDir","normal");
        xlabel('u'); ylabel('v'); title('Input Image');
    subplot(1,6,2); imshow(grid_img'); set(gca,"YDir","normal");
        title('Start Coordinate Grid');
    subplot(1,6,3); imshow(search_img');
        set(gca,"YDir","normal"); title('Raw Detection Result');
    subplot(1,6,4); imshow(search_finalImg'); set(gca,"YDir","normal");
        title('Post-Processing Result')
    subplot(1,6,5); imshow(permute(labels_c,[2 1 3]));
        set(gca,"YDir","normal"); title('Labeled ROIs');
    subplot(1,6,6); imshow(found_branches'); set(gca,"YDir","normal");
        title({'Found Branch ROIs','Home ROI Removed'});
%}
end