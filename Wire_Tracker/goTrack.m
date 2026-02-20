%% main tracking logic from start P with vector N
% Jade Lariviere | last modified Nov. 18, 2025

function [Run] = goTrack(Wire,P,N,Opts,Param,Loop)
% this function orchestrates the main tracking logic using the information
% of the "Wire" (a presumably tube-like shape, but can be anything), a
% start point and a start direction. The tracking behavior and duration is 
% adjusted based on the Parameters and Loop structure settings. The Run 
% structure is returned, containing relevant information including the raw 
% and smooth tracked centerline and possible branches found for further
% tracking iterations. ====================================================

% initialize loop & certain variables
tag = Opts.gifDir + string(datetime,'yyyyMMMdd_HH-mm-ss'); % recording name
windowState = Opts.figView{2}; % initial windowState of figure in loop
i = 1; % initialize loop counter
num_flips = 0; % initialize flips counter
dims = Wire.dims; % wire volume dimensions
N = N/vecnorm(N); % ensure N is normalized

Centerline = nan([Loop.MaxTrackIt+1, 3]); % tracked coordinates
    Centerline(1,:) = P;
Direction = nan([Loop.MaxTrackIt+1, 3]); % tracked direction vectors
    Direction(1,:) = N;
Run.Branch.area = zeros(size(Wire.volume)); % store detected branches
Run.Branch.dirs = cell(size(Wire.volume)); % store branch directions

% now, GO TRACK! ==========================================================
while (i < Loop.MaxTrackIt && num_flips < Loop.MaxFlipFlags)
P_old = P; N_old = N;
% determine furthest direction, then step ---------------------------------
[Cam.Data,arePts] = getProjection(P,N,Wire,'Ahead',Param.AheadPlaneTol);
if ~arePts % seemingly found self peering outside of volume 
    warning("goTrack() stopped, had empty view."); break;
end
[Cam.Image] = makeImage(Cam.Data,Param.ImageDims,'Ortho');
[N_new,step_size,Cam.FinImage] = getOrthoInfo(Cam.Data,Cam.Image,'Absolute');
    N = N_old*Param.Momentum + N_new*(1 - Param.Momentum);
        N = N/vecnorm(N);
    P = P_old + N*(step_size*Param.Speed); % ------------------------------

% settle into center of new cross section % -------------------------------
[Slice.Data,arePts] = getProjection(P,N,Wire,'On',Param.OnPlaneTol);
if ~arePts % seemingly found self outside of volume 
    warning("goTrack() stopped after step."); break;
end
[Slice.Image] = makeImage(Slice.Data,Param.ImageDims,'Slice');
[P_new,Slice.FinImage] = getVolumeInfo(Slice.Data,Slice.Image,'Local');
    P = P_new; % ----------------------------------------------------------

% store possible branching + directions from slice ------------------------
if Opts.checkBranches
    [Run.Branch.area,Run.Branch.dirs] = findBranch(Run.Branch.area,Slice, ...
        Run.Branch.dirs,i);
end

% plotting for visualization ----------------------------------------------
if Opts.viewTrack
    figure(2); set(gcf,"WindowState",windowState);
    subplot(2,2,1:2:3);
        scatter3(Wire.verts(1:2:end,1),Wire.verts(1:2:end,2), ...
            Wire.verts(1:2:end,3),'k.','LineWidth',0.1);
        hold on; scatter3(P(1),P(2),P(3),'r','filled','LineWidth',3);
        quiver3(P(1),P(2),P(3),N(1)*5,N(2)*5,N(3)*5,'off','r','LineWidth',2);
        hold off; view(3); xlabel('x'); ylabel('y'); zlabel('z');
        title("Tracking Process"); view(Opts.figView{1}); axis equal;
    subplot(2,2,2);
        imshow(Cam.FinImage.ROI',Cam.FinImage.range);
        set(gca,"YDir","normal"); xlabel('u'); ylabel('v'); colormap("pink");
        title("Camera");
    subplot(2,2,4);
        imshow(Slice.FinImage.ROI',Slice.FinImage.range);
        set(gca,"YDir","normal"); xlabel('u'); ylabel('v'); drawnow;
        title("Slice");
    if Opts.makeRunGif % produce GIF of wire tracking ---------------------
        exportgraphics(figure(2),sprintf("%s.gif",tag),'Append',true);
    end
end

% save point & direction, update loop counter -----------------------------
Centerline(i+1,:) = P;      Direction(i+1,:) = N;
i = i + 1;
% -------------------------------------------------------------------------
% determine approximately how close tracker is to wall
dist_close = median(Cam.FinImage.ROI(Cam.FinImage.ROI~=0),'all');
% terminate loop early on some conditions
P_round = [min(ceil(P(1)),dims(1)), min(ceil(P(2)),dims(2)), ...
    min(ceil(P(3)),dims(3))];
if Wire.volume(P_round(1),P_round(2),P_round(3)) == 0
    fprintf('Tracking exited wire volume. Aborting...\n'); break;
elseif norm(dist_close) < Param.ExitTol
    fprintf('Approaching end (surface wall). Aborting...\n'); break;
end
% check for tracker flip --------------------------------------------------
if (i ~= 1) && (rem(i,Param.Memory)==1)
%if i > Param.Memory
    % pass last points within memory
    [hasFlipped, P_new, N_new] = flipCheck(Centerline(i-Param.Memory:i,:));
    if hasFlipped
        num_flips = num_flips + 1; % add to counter
        fprintf('Flip %g of %g flagged (i = %g).\n',num_flips,Loop.MaxFlipFlags,i);
        % erase previous points & direction based on Param.Memory
        Centerline(i-Param.Memory:i,:) = nan;
        Direction(i-Param.Memory:i,:) = nan;
        P = P_new; N = N_new; % assign (revert) new P and N
    end
end
windowState= get(gcf,"WindowState"); % update state if user changed it
end

% end tracking, clean up data =============================================
fprintf('Tracking iteration %g of %g reached.\n',i,Loop.MaxTrackIt);
if Opts.viewTrack, close gcf; end

% prune centerline list of missing data
Run.Points.raw = rmmissing(Centerline);
Run.Vectors = rmmissing(Direction);
% average adjacent points for smoothing
pts_next = [Run.Points.raw(2:end,:); Run.Points.raw(end,:)];
pts_before = [Run.Points.raw(1,:); Run.Points.raw(1:end-1,:)];
Run.Points.smooth = (Run.Points.raw + pts_before + pts_next)./3;
% store Run info
Run.numIterations = i;
Run.numPts = size(Run.Points.smooth,1);
Run.numFlips = num_flips; % num iterations =/= num final points

% clean up branching detection volume -------------------------------------
if Opts.checkBranches
    se = strel("sphere",2); % arbitrary structural element
    Run.Branch.ROIs = imdilate(Run.Branch.area,se); % lump nearby ROIs
        Run.Branch.ROIs = Run.Branch.ROIs.*(Wire.volume~=0); % crop to vol
        [Run.Branch.ROIs,num_ROIs] = bwlabeln(Run.Branch.ROIs,18);
    % return Tracker structure: branch ROI centroids
    branch_centroids = regionprops(Run.Branch.ROIs,'Centroid');
    branch_Ps = cell2mat(struct2cell(branch_centroids)');
    if num_ROIs == 0
        fprintf('----- goTrack did not detect branching this run. -----\n');
        Run.Branch.nextPs = []; Run.Branch.nextNs = [];
        Run.Branch.numBranch = 0; return;
    end
    % X & Y are swapped to preserve column organization in rest of code
    Run.Branch.nextPs = [branch_Ps(:,2) branch_Ps(:,1) branch_Ps(:,3)];

    % return Tracker structure: branch ROI average vector
    Run.Branch.nextNs = cell(num_ROIs,1);
    for i = 1:num_ROIs
        vector_idx = find(Run.Branch.ROIs==i);
        branch_vectors = Run.Branch.dirs(vector_idx);
            branch_vectors = branch_vectors(~cellfun('isempty',branch_vectors));
        if any(~cellfun('isempty',branch_vectors)) % if a vector exists
            Run.Branch.nextNs{i} = mean(vertcat(branch_vectors{:}),1);
        else % else, create a random vector in case
            Run.Branch.nextNs{i} = [0 0 1];
        end
    end
    Run.Branch.nextNs = cell2mat(Run.Branch.nextNs);
    Run.Branch.numBranch = num_ROIs; % store info of num of next pts
else, Run.Branch.numBranch = 0;
end
end