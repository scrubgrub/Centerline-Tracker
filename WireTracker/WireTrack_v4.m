%% WireTracker v4
% coded by Jade Lariviere
% last edited: Mar. 23, 2025

%% PRE-PROCESSING =========================================================
clc; clear; close all;
% EDITABLE PARAMETERS =====================================================
Param.OnPlaneTol    = 0.5; % volume projection; larger = more forgiving
Param.AheadPlaneTol = 0.01; % ortho camera; larger = more 'farsighted'
Param.ExitTol       = 1.7; % abort before hitting wall; larger = farther
Param.SurfTol       = 1e-2; % isosurface value; smaller = expands surface
Param.Momentum      = 0.3; % how much of previous N weighted in N_new
Param.Speed         = 0.8; % multiplier for speed; higher = faster jumps
Param.ImageDims     = [70 70]; % returned image dimensions
Param.erosion       = 2; % structural element size for erosion
Param.smoothing     = 5; % amount volume is blurred for weighted volume
Param.pad           = 10; % number of pixels to pad volume for processing

Loop.MaxSearchIt    = 80; % max number of search iterations
Loop.MaxTrackIt     = 171; % max number of track iterations (< 171!!)

% =========================================================================
% import file info % do preliminary processing...
[I.name,I.path] = uigetfile({'*.raw;*.nii.gz;*.nii','Structural Data'}, ...
    'Select file with wire segmentation data.');
    [~,~,I.ext] = fileparts(I.name); % get input file info
if strcmp(I.ext,'.raw') % for .raw files
    Wire.raw = double(multibandread(strcat(I.path,'\',I.name), ...
    [224 224 224],'*uint8',0,'bsq','ieee-be')); % fill file info by hand
        Wire.raw = Wire.raw./max(Wire.raw,[],'all'); % 0 - 1 best results
else % for .nii.gz files
    Wire.raw = double(niftiread(strcat(I.path,'\',I.name)));
end

I.dims = size(Wire.raw); Wire.dims = I.dims; % file dimensions
vol_padded = padarray(Wire.raw,[Param.pad Param.pad Param.pad]); % pad vol
[x,y,z] = ndgrid(1:I.dims(1),1:I.dims(2),1:I.dims(3)); % create grid

% creating a weighted volume; blurring = center > edges
se = strel('sphere',Param.erosion); % structural element for erosion
wire_eroded = imerode(vol_padded,se); % eroded image
Wire.volume = smooth3(wire_eroded,'gaussian',5,Param.smoothing);
    Wire.volume = Wire.volume(Param.pad+1:end-Param.pad, ...
              Param.pad+1:end-Param.pad,Param.pad+1:end-Param.pad);

% creating wire shell using isosurface() -> convert to 3D volume
[Wire.faces,Wire.all_verts] = isosurface(x,y,z,Wire.volume,Param.SurfTol);
Wire.verts = unique(round(Wire.all_verts),'rows'); % preserve unique verts
idx_surf = sub2ind(I.dims,Wire.verts(:,1),Wire.verts(:,2),Wire.verts(:,3));
Wire.surface = zeros(I.dims);
    Wire.surface(idx_surf) = 1;

% =========================================================================
% create a patch from isosurface for UI stuff
Wire.Patch.object = makePatch(Wire,[0.5 0.6 0.9]); set(gcf,'visible','off');

% request to proceed automatically or not
answer = questdlg('Run AUTO END SEARCH, or INPUT START POINT & DIRECTION?',...
    'Start Determination','Automatic','Manual Entry','Automatic');
switch answer
case 'Automatic' % --------------------------------------------------------
    % give algorithm start coordinate for tracking = max weighted volume pt
    [~,idx_Pmax] = max(Wire.volume,[],'all');
    [P_x,P_y,P_z] = ind2sub(Wire.dims,idx_Pmax);
        P = [P_x, P_y, P_z];
    % give algorithm start direction via vector between max and 2nd max
    [~,idx_Pmax2] = max(Wire.volume(Wire.volume<max(Wire.volume)));
    [P_x,P_y,P_z] = ind2sub(Wire.dims,idx_Pmax2);
        N = [P_x, P_y, P_z] - P;
        N = N/norm(N);
case 'Manual Entry' % -----------------------------------------------------
        f = figure; f_ax = axes; % new figure for display to user
        set(gcf,'visible','off');
        hold on; copyobj(Wire.Patch.object,f_ax);
        scatter3(Wire.verts(:,1),Wire.verts(:,2),Wire.verts(:,3),'k.');
        hold off;
    answer = dialogWireStart(f); % open dialog box
    P = answer(1,:); N = answer(2,:); % set starting points & direction!
        N = N/norm(N); % normalize N
end % ---------------------------------------------------------------------

clear se x y z *_*
% continue to next section to actually run the wire!
%% LOOP 1: AUTOMATIC WIRE END SEARCH ======================================
% find wire end using one starting location in wire
if strcmp(answer,'Automatic') % ONLY IF 'AUTOMATIC' WAS CHOSEN!
fprintf('--- STARTING WIRE END SEARCH... ---\n');
Loop.num = 1; % initialize loop counter
Wire.searchPts = nan([Loop.MaxSearchIt+1 3]); % store search coordinate pts
    Wire.searchPts(1,:) = P;
Wire.searchDirs = nan([Loop.MaxSearchIt+1 3]); % store search directions
    Wire.searchDirs(1,:) = N;

while (Loop.num < Loop.MaxSearchIt)
P_old = P; N_old = N;
% determine furthest direction, then step ---------------------------------
Ortho.Data = getProjection(P,N,Wire,'Ahead',Param.AheadPlaneTol);
[Ortho.Image] = makeImage(Ortho.Data,Param.ImageDims,'Ortho');
[N_new,step_size,Ortho.FinImage] = getOrthoInfo(Ortho.Data,Ortho.Image,'Absolute');
    N = N*Param.Momentum + N_new*(1 - Param.Momentum);
    P = P + N_new*(step_size*Param.Speed); % -------------------------------

% settle into center of new cross section % -------------------------------
Slice.Data = getProjection(P,N,Wire,'On',Param.OnPlaneTol);
[Slice.Image] = makeImage(Slice.Data,Param.ImageDims,'Slice');
[P_new,Slice.FinImage] = getVolumeInfo(Slice.Data,Slice.Image,'Local');
    P = P_new; % ----------------------------------------------------------

% plotting for visualization ----------------------------------------------
figure(2); 
    scatter3(Wire.verts(:,1),Wire.verts(:,2),Wire.verts(:,3),'.');
    hold on; scatter3(P(1),P(2),P(3),'r','filled','LineWidth',3);
    quiver3(P(1),P(2),P(3),N(1)*5,N(2)*5,N(3)*5,'off','r','LineWidth',2);
    hold off; view(3);
    xlabel('x'); ylabel('y'); zlabel('z');
drawnow; %pause(0.1);

% save point & direction history, update loop counter ---------------------
Wire.searchPts(Loop.num+1,:) = P_new;
Wire.searchDirs(Loop.num+1,:) = N_new;
Loop.num = Loop.num + 1;
% -------------------------------------------------------------------------
% find shortest distance of point on wire edge to P
dist_close = min(Ortho.Data.distance,[],'all');
% terminate loop early on some conditions
if Wire.raw(round(P(1)),round(P(2)),round(P(3))) == 0
    fprintf('Tracking exited wire volume. Aborting...\n'); break;
elseif norm(dist_close) < Param.ExitTol
    fprintf('Approaching end (surface wall). Aborting...\n'); break;
end
end % ---------------------------------------------------------------------
fprintf('Search iteration %g of %g reached.\n',Loop.num,Loop.MaxSearchIt);
close gcf;

% =========================================================================
% remove missing entries
Wire.searchPts = rmmissing(Wire.searchPts);
Wire.searchDirs = rmmissing(Wire.searchDirs);

% use last 3 P and N to get approximate starting location + direction
P = mean(Wire.searchPts(end-2:end,:),1);
N = -mean(Wire.searchDirs(end-2:end,:),1); % invert to go opposite dir
end

% continue to next section to track entirety of wire!
%% LOOP 2: START FROM WIRE END AND GO TO OTHER END ========================
% begin main wire tracking
fprintf('--- STARTING WIRE TRACKING... ---\n');
Loop.num = 1; % initialize loop counter
Centroid.raw = nan([Loop.MaxTrackIt+1 3]); % store search coordinate pts
    Centroid.raw(1,:) = P;

while (Loop.num < Loop.MaxTrackIt)
P_old = P; N_old = N;
% determine furthest direction, then step ---------------------------------
Ortho.Data = getProjection(P,N,Wire,'Ahead',Param.AheadPlaneTol);
[Ortho.Image] = makeImage(Ortho.Data,Param.ImageDims,'Ortho');
[N_new,step_size,Ortho.FinImage] = getOrthoInfo(Ortho.Data,Ortho.Image,'Absolute');
    N = N*Param.Momentum + N_new*(1 - Param.Momentum);
    P = P + N_new*(step_size*Param.Speed); % -------------------------------

% settle into center of new cross section % -------------------------------
Slice.Data = getProjection(P,N,Wire,'On',Param.OnPlaneTol);
[Slice.Image] = makeImage(Slice.Data,Param.ImageDims,'Slice');
[P_new,Slice.FinImage] = getVolumeInfo(Slice.Data,Slice.Image,'Local');
    P = P_new; % ----------------------------------------------------------

% plotting for visualization ----------------------------------------------
figure(2);
subplot(2,2,1:2:3);
    scatter3(Wire.verts(:,1),Wire.verts(:,2),Wire.verts(:,3),'.');
    hold on; scatter3(P(1),P(2),P(3),'r','filled','LineWidth',3);
    quiver3(P(1),P(2),P(3),N(1)*5,N(2)*5,N(3)*5,'off','r','LineWidth',2);
    hold off; view(3);
    xlabel('x'); ylabel('y'); zlabel('z');
subplot(2,2,2);
    imshow(Ortho.FinImage.ROI',[0 50]); set(gca,"YDir","normal");
    xlabel('u'); ylabel('v'); colormap("pink");
subplot(2,2,4);
    imshow(Slice.FinImage.ROI',[0 0.5]); set(gca,"YDir","normal");
drawnow; %pause(0.1);
%exportgraphics(figure(1),'sub004OZtrack.gif','Append',true);

% save centroid, update loop counter --------------------------------------
Centroid.raw(Loop.num+1,:) = P_new;
Loop.num = Loop.num + 1;
% -------------------------------------------------------------------------
% find shortest distance of point on wire edge to P
dist_close = min(Ortho.Data.distance,[],'all');
% terminate loop early on some conditions
if Wire.raw(round(P(1)),round(P(2)),round(P(3))) == 0
    fprintf('Tracking exited wire volume. Aborting...\n'); break;
elseif norm(dist_close) < Param.ExitTol
    fprintf('Approaching end (surface wall). Aborting...\n'); break;
end
end % ---------------------------------------------------------------------
fprintf('Tracking iteration %g of %g reached.\n',Loop.num,Loop.MaxTrackIt);
close gcf;

% continue to next section to perform post-processing of data!
%% POST PROCESSING ========================================================
% clean up wire centroids (prune list & average points)
Centroid.raw = rmmissing(Centroid.raw);
pts_next = [Centroid.raw(2:end,:); Centroid.raw(end,:)]; % averaging
pts_before = [Centroid.raw(1,:); Centroid.raw(1:end-1,:)];
    Centroid.smooth = (Centroid.raw + pts_before + pts_next)./3;

f = figure; f_ax = axes; set(gcf,'visible','off'); % figure for dialog box -----------------
    hold on; copyobj(Wire.Patch.object,f_ax);
    plot3(Centroid.smooth(:,1),Centroid.smooth(:,2),Centroid.smooth(:,3),'.-');
    scatter3(Centroid.smooth(1,1),Centroid.smooth(1,2),Centroid.smooth(1,3),'cyan','filled');
    scatter3(Centroid.smooth(end,1),Centroid.smooth(end,2),Centroid.smooth(end,3),'red','filled');
    hold off;
answer = dialogTrackingCheck(f); % open dialog for final check
if answer; Centroid.smooth = flipud(Centroid.smooth); end % flip if yes

clear *_*
% continue to next section to export your data!
%% SAVE & EXPORT ==========================================================
% save .mat file of centroid array!
O.path = uigetdir(I.path,'Select folder to save wire centroids.');
    O.path = strrep(O.path,'\','/'); % change slashes

% get user input to name file, generate name, and save --------------------
answer = inputdlg({'Subject ID','Volume Identifier'},'Name .mat file.',...
    [1 15; 1 15]);
O.name = sprintf('/%s_%s.mat',answer{1},answer{2});

save(strcat(O.path,O.name),'-struct','Centroid','smooth');
fprintf('Saving %s to %s...\n',O.name,O.path); % --------------------------

clear f answer *_*
% and we're done!
