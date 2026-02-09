%% helper function to acquire specific slice data for further use
% Jade Lariviere | last modified Dec. 7, 2025

function [sliceData,isPoints] = getProjection(P,N,volumeData,ProjectionType,planeTol)
% function takes argument volumeData which contains info about the volume 
% size and values and returns structure sliceData with: points [u v] within 
% the planeTol threshold, their off-plane separation, and other useful 
% slice information. requires a point P [x,y,z] and normal vector N [x,y,z] 
% to define the plane. the input argument ProjectionType determines the 
% comparison method the function will perform. ============================

dims = volumeData.dims;
N_hat = N/norm(N); % normalize vector N in case
frac = 3; % optimization; only look at 1/n of the total volume for voxels

switch ProjectionType % select which data type to use
    case 'On'
        % sample only a portion of total volume voxels in data
        bound_x = ceil(max(P(1)-dims(1)/frac,1):min(P(1)+dims(1)/frac,dims(1)));
        bound_y = ceil(max(P(2)-dims(2)/frac,1):min(P(2)+dims(2)/frac,dims(2)));
        bound_z = ceil(max(P(3)-dims(3)/frac,1):min(P(3)+dims(3)/frac,dims(3)));
        data = volumeData.volume;
        sample = data(bound_x,bound_y,bound_z);
        sample_bound = [bound_x(1), bound_y(1), bound_z(1)] - 1; % get corner
        indices = find(sample); % find nonzeros
    case 'Ahead'
        % no subsampling from data, no offset
        data = volumeData.surface; sample = data;
        sample_bound = [0, 0, 0];
        indices = volumeData.idxSurf; % find nonzeros
    otherwise
        error('Input argument %s for ProjectionType not recognized.', ...
            ProjectionType);
end

% nonzero coordinates from sample indices =================================
[ptX, ptY, ptZ] = ind2sub(size(sample),indices); % find nonzeros
ptXYZ = [ptX+sample_bound(1) ptY+sample_bound(2) ptZ+sample_bound(3)]; % reformat
    n = repmat(N_hat,[length(ptXYZ),1]); % for vector math

% plane logic =============================================================
switch ProjectionType
    case 'On' % select points ON the plane
        isPtGood = abs(dot(P-ptXYZ,n,2)) < planeTol; % dot(ortho vecs) = 0
    case 'Ahead' % select points AHEAD of plane
        isPtGood = dot(P-ptXYZ,n,2) < planeTol; % dot(ortho vecs) = 0
end

if any(isPtGood) % check if ANY points exist on plane & flag
    isPoints = true;
else
    warning("getProjection() found no points %s plane.\n",ProjectionType);
    isPoints = false; sliceData = []; return;
end

% return list of point coordinates ========================================
good_pts = ptXYZ(isPtGood,:); % point coordinates
good_vals = data(sub2ind(dims,good_pts(:,1),good_pts(:,2),good_pts(:,3)));

% project 3D points onto plane's coord space ==============================
% origin of point P with ortho vectors of N as the axes
orthoNs = null(N_hat);
axisU = orthoNs(:,1)'; axisV = orthoNs(:,2)';
numPts = size(good_pts,1); % dupe N, U and V for vectorized math
    n = repmat(N_hat,[numPts,1]);
    aU = repmat(axisU,[numPts,1]);
    aV = repmat(axisV,[numPts,1]);

% transform [x,y,z] into new 2D coordinate space [u,v]
u_raw = dot(aU,good_pts-P,2); v_raw = dot(aV,good_pts-P,2);
s_raw = dot(n,good_pts-P,2); % out of plane separation
    s_raw(s_raw < 0) = 0; % set any negative values to zero
dist_raw = vecnorm(good_pts-P,2,2); % distance from point to P

[s_sorted, idx_sorted] = sort(s_raw,'descend'); % sort furthest -> closest
u_sorted = u_raw(idx_sorted,:); v_sorted = v_raw(idx_sorted,:);
dist_sorted = dist_raw(idx_sorted);
values_sorted = good_vals(idx_sorted);

% return sliceData structure ==============================================
sliceData.pointsUV = [u_sorted v_sorted];
    sliceData.pointsXYZ = good_pts(idx_sorted);
    sliceData.separation = s_sorted;
    sliceData.distance = dist_sorted;
    sliceData.intensity = values_sorted;
    sliceData.P = P;            sliceData.N = N_hat;
    sliceData.axes.U = axisU;   sliceData.axes.V = axisV;

% plots for debugging =====================================================
%{
figure(10);
subplot(1,2,1);
    scatter3(ptXYZ(:,1),ptXYZ(:,2),ptXYZ(:,3),'k.'); hold on;
    scatter3(good_pts(:,1),good_pts(:,2),good_pts(:,3),'filled');
    scatter3(P(1),P(2),P(3),'blue','filled');
    hold off; view(3); xlabel('x'); ylabel ('y'); zlabel('z');
    title('Found Points');
subplot(1,2,2); scatter(u_raw,v_raw,'filled'); xlabel('u'); ylabel('v');
    title('Projected Points in UV Space');
%}
end
