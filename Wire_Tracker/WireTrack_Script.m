%% WireTracker v2.0
% coded by Jade Lariviere
% last edited: Dec. 8, 2025
    
%% PRE-PROCESSING =========================================================
clc; clear; close all;
% OPTIONS & PARAMETERS ====================================================
Opts.viewTrack      = false; % (slows!) toggle to view tracking process;
Opts.checkBranches  = true; % toggle to permit branch detection/tracking
Opts.doInterp       = true; % interpolate branch paths to main centerline
Opts.makeRunGif     = false; % (slows!) if viewTrack = true, record 1st run

Param.OnPlaneTol    = 0.5; % cross section; larger = more forgiving
Param.AheadPlaneTol = 0.05; % camera proximity thresh; larger = 'farsighted'
Param.ExitTol       = 1.45; % end tracking before a wall; larger = farther  
Param.Momentum      = 0.10; % how much of previous N weighted in N_new
Param.Speed         = 0.30; % multiplier for speed; higher = faster jumps
Param.Memory        = 6; % n previous points passed to flipCheck (>= 3)
Param.ImageDims     = [65 65]; % returned image dimensions
Param.smoothing     = 0; % volume smoothing parameter (SD of guassian blur)
Param.pad           = 5; % number pixels to pad volume for smoothing
Param.PruneLength   = 3; % number points or fewer in a branch to prune

Loop.MaxTrackIt     = 250; % max number of track iterations
Loop.MaxFlipFlags   = 3; % max number of allotted flips per tracking run

% =========================================================================
[I, Wire] = openFile;

% extract Wire volume + surface -------------------------------------------
if Param.smoothing == 0, Wire.volume = Wire.raw; % wire interior 
    
else
    smooth_volume = imgaussfilt3(Wire.raw,Param.smoothing);
    Wire.volume = smooth_volume; % wire interior 
end
% wire surface = voxels just outside of vol perimeter + boundary voxels
Wire.surface = imdilate(Wire.volume,strel("cuboid",[3 3 3]));
    Wire.surface = Wire.surface.*(~Wire.volume) + Wire.volume.*I.bounds;
Wire.idxSurf = find(Wire.surface); % preallocate known nonzero voxels
fprintf('Identified structure volume and surface.\n');

% create visualization of volume geometry using isosurface() for UI
vol_padded = padarray(Wire.volume,[Param.pad Param.pad Param.pad]); % pad vol
    pad_dims = size(vol_padded);
[x,y,z] = ndgrid(1:pad_dims(1),1:pad_dims(2),1:pad_dims(3)); % create grid
[Wire.faces,Wire.all_verts] = isosurface(x,y,z,vol_padded,0);
    Wire.all_verts = Wire.all_verts - Param.pad; % remove coordinate offset
Wire.verts = unique(round(Wire.all_verts),'rows'); % keep unique verts

% =========================================================================
% create a patch from isosurface for UI stuff
Wire.Patch.object = makePatch(Wire,[0.5 0.6 0.9]); set(gcf,'visible','off');
fprintf('Isosurface object created.\n')

% request to proceed automatically or not
answer = questdlg('Run AUTO END SEARCH, or INPUT START POINT & DIRECTION?',...
    'Start Determination','Automatic','Manual Entry','Automatic');
switch answer
case 'Automatic' % --------------------------------------------------------
    % create approximations of end points with bwskel + bwmorph3
    wire_skel = bwskel(logical(Wire.volume),"MinBranchLength",10);
    wire_ends = bwmorph3(wire_skel,"endpoints");
    [P_x, P_y, P_z] = ind2sub(size(wire_ends),find(wire_ends));
    f = figure("Visible","off"); f_ax = axes; % UI fig
        hold on; copyobj(Wire.Patch.object,f_ax);
        scatter3(P_x,P_y,P_z,'ro','filled');
        hold off;
    [P_labels] = num2cell(1:size(P_x,1));
        label_offset = 1;
    text(P_x,P_y,P_z+label_offset,P_labels,"FontSize",9,'FontWeight','bold');
    [answer,Opts.figView] = dialogWireStart(f,[],[P_x,P_y,P_z]);
case 'Manual Entry' % -----------------------------------------------------
        f = figure; f_ax = axes; % new figure for display to user
        set(gcf,'visible','off');
        hold on; copyobj(Wire.Patch.object,f_ax);
        scatter3(Wire.verts(:,1),Wire.verts(:,2),Wire.verts(:,3),'k.');
        hold off;
    [answer,Opts.figView] = dialogWireStart(f,[]); % open dialog box
end % ---------------------------------------------------------------------
P = answer(1,:); N = answer(2,:); % set starting points & direction!

% ensure dir for storing media if not already -----------------------------
Opts.gifDir = fileparts(mfilename('fullpath')) + "\recordings\";
if Opts.makeRunGif && ~exist(Opts.gifDir,"dir"), mkdir(Opts.gifDir); end

clear se x y z *_* answer
% continue to next section(s) to actually run the wire!
%% PHASE 1: START TRACKING UNTIL REACHING AN END ==========================
% begin main wire tracking, initialize relevant variables
fprintf('=============== STARTING WIRE TRACKING... ===============\n');
%Opts.checkBranches = setting_branch;
tic;    Run = goTrack(Wire,P,N,Opts,Param,Loop);

% continue to next section to perform additional processing of data!
%% PHASE 2: BRANCH EXPLORATION ITERATIVE LOOPING ==========================
% create struct array for all Tracker calls
field_names = fieldnames(Run);
num_fields = size(field_names,1); num_branch = Run.Branch.numBranch;
% create struct to store all tracking iterations
    Tracker = cell2struct(cell([num_fields,num_branch+1]),field_names);
        Tracker(1) = Run;

if num_branch ~= 0 % tracking loop ----------------------------------------
    % num_Ns should be equal to num_Ps
    next_P = Run.Branch.nextPs; next_N = Run.Branch.nextNs;
        setting_checkBranch = Opts.checkBranches;
        Opts.checkBranches = false; Opts.viewTrack = false;
    for i = 1:num_branch % loop through all starting P and N!
        fprintf('------- Tracking branch %d of %d. -------\n',i,num_branch);
        Tracker(i+1) = goTrack(Wire,next_P(i,:),next_N(i,:),Opts,Param,Loop);
    end
    Opts.checkBranches = setting_checkBranch; % reset to setting
end

clear *_*
% continue to the next section to finalize results
%% POST PROCESSING ========================================================
% prune short branches & sort order before finalizing and export
if Opts.checkBranches
    branch_sizes = arrayfun(@(x) size(x.raw,1),[Tracker.Points]);
    short_branches = branch_sizes <= Param.PruneLength;
        short_branches(1) = 0; % first branch prevented from being pruned!
    Tracker(short_branches) = [];
    sort_near = arrayfun(@(x) norm(x.Points.smooth(1,:)-Tracker(1).Points.smooth(1,:)), ...
    Tracker,'UniformOutput',false);
    [~,sort_idx] = sort(cell2mat(sort_near));
    Tracker = Tracker(sort_idx);
    fprintf(2,"*** Number of branches found: %d\n" + ...
        "*** Branches pruned: %d\n*** Final branch number: %d\n", ...
        Tracker(1).Branch.numBranch,sum(short_branches),size(Tracker,1)-1);
end

% finalize the tracking paths ---------------------------------------------
Tracker(1).Points.final = Tracker(1).Points.smooth; % first curve
if Opts.checkBranches && Opts.doInterp % remaining tracks interpolated
    for i = 2:length(Tracker)
        Tracker(i).Points.final = interpBranch(Tracker(i).Points.smooth,...
            Tracker(1).Points.final,Param);
    end
else % pass remaining smoothed tracks as final
    for i = 2:length(Tracker)
        Tracker(i).Points.final = Tracker(i).Points.smooth;
    end
end;    toc;

% dialog box figure -------------------------------------------------------
f = figure("WindowState","maximized","Visible","off"); f_ax = axes; % UI fig
    hold on; copyobj(Wire.Patch.object,f_ax);
    plot3(Tracker(1).Points.final(:,1),Tracker(1).Points.final(:,2), ...
            Tracker(1).Points.final(:,3),'k.-','LineWidth',1.5);
    arrayfun(@(x) scatter3(x.Points.smooth(1,1),x.Points.smooth(1,2), ...
        x.Points.smooth(1,3), ...
        'cyan','filled'), Tracker); % start points
    arrayfun(@(x) scatter3(x.Points.smooth(end,1),x.Points.smooth(end,2), ...
        x.Points.smooth(end,3), ...
        'red','filled'), Tracker); % end points
    for i = 2:size(Tracker,1)
    plot3(Tracker(i).Points.final(:,1),Tracker(i).Points.final(:,2), ...
        Tracker(i).Points.final(:,3),'.-');
    text(Tracker(i).Points.smooth(1,1),Tracker(i).Points.smooth(1,2), ...
        Tracker(i).Points.smooth(1,3), ...
        sprintf('Brc %g',i-1),"FontSize",9); % label
    end
    hold off;
% final view dialog, flip if yes
track_length = [Tracker.numPts];
answer = dialogTrackCheck(f,Opts.figView,track_length);
if answer; Tracker(1).Points.smooth = flipud(Tracker(1).Points.smooth); end

clear *_* se answer
% continue to next section to export your data!
%% SAVE & EXPORT ==========================================================
% save .mat file of centroid array!
O.path = uigetdir(I.path,'Select folder to save wire centroids.');
    O.path = strrep(O.path,'\','/'); % change slashes

% get user input to name file, generate name, and save --------------------
answer = inputdlg({'Subject ID','Volume Identifier'},'Name .mat file.',...
    [1 15; 1 15]);
O.name = sprintf('/%s_%s.mat',answer{1},answer{2});
O.centerlines = cell2struct(arrayfun(@(x) x.Points.final,Tracker,'UniformOutput',false), ...
    'Coordinates',2);

save(strcat(O.path,O.name),'-struct','O','centerlines');
fprintf(2,'Saved %s to %s !\n',O.name,O.path); % --------------------------

clear f answer *_*
% and we're done!
