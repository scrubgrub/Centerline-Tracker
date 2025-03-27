%% helper function to acquire specific slice data for further use
% Jade Lariviere | last modified Mar. 5, 2025

function [sliceData] = getProjection(P,N,volumeData,ProjectionType,planeTol)
% function takes input structure volumeData which contains info about the
% volume size and values (size m x n x p) and returns a structure
% sliceData with the projection of points [u v], separation, intensity, and
% slice information. requires a point P [x,y,z] and normal vector N [x,y,z] 
% to define the plane. the input argument ProjectionType determines what 
% type of comparison the function will perform. ===========================

dims = volumeData.dims;
N_unit = N/norm(N); % normalize vector N
switch ProjectionType % select which data type to use
    case 'On'; data = volumeData.volume; 
    case 'Ahead'; data = volumeData.surface;
    otherwise
        error('Input argument %s for ProjectionType not recognized.', ...
            ProjectionType);
end

% grab nonzero coordinates from 3D volume =================================
[ptX, ptY, ptZ] = ind2sub(dims,find(data)); % find nonzeros
ptXYZ = [ptX ptY ptZ]; % coords of wire points
    n = repmat(N_unit,[length(ptXYZ),1]); % for vector math

% plane logic =============================================================
switch ProjectionType
    case 'On' % select points ON the plane
        isPtGood = abs(dot(P-ptXYZ,n,2)) < planeTol; % dot(ortho vecs) = 0
    case 'Ahead' % select points AHEAD of plane
        isPtGood = dot(P-ptXYZ,n,2) < planeTol; % dot(ortho vecs) = 0
end

% return list of point coordinates ========================================
ptsGood = ptXYZ(isPtGood,:); % point coordinates
ptsValues = data(sub2ind(dims,ptsGood(:,1),ptsGood(:,2),ptsGood(:,3)));

% project 3D points onto plane's coord space ==============================
% origin of point P with ortho vectors of N as the axes
orthoNs = null(N_unit);
axisU = orthoNs(:,1)'; axisV = orthoNs(:,2)';
numPts = length(ptsGood); % dupe N, U and V for vectorized math
    n = repmat(N_unit,[numPts,1]);
    aU = repmat(axisU,[numPts,1]);
    aV = repmat(axisV,[numPts,1]);

% transform [x,y,z] into new 2D coordinate space [u,v]
u_raw = dot(aU,ptsGood-P,2); v_raw = dot(aV,ptsGood-P,2);
s_raw = dot(n,ptsGood-P,2); % out of plane separation
dist_raw = vecnorm(ptsGood-P,2,2); % distance from point to P

[s_sorted, idx_sorted] = sort(s_raw,'descend'); % sort furthest -> closest
u_sorted = u_raw(idx_sorted,:); v_sorted = v_raw(idx_sorted,:);
dist_sorted = dist_raw(idx_sorted);
values_sorted = ptsValues(idx_sorted);

% return sliceData structure ==============================================
sliceData.pointsUV = [u_sorted v_sorted];
    sliceData.pointsXYZ = ptXYZ;
    sliceData.separation = s_sorted;
    sliceData.distance = dist_sorted;
    sliceData.intensity = values_sorted;
    sliceData.P = P;            sliceData.N = N;
    sliceData.axes.U = axisU;   sliceData.axes.V = axisV;

% plots for debugging =====================================================
%{
figure;
subplot(1,2,1); hold on;
    scatter3(ptX,ptY,ptZ,'k.'); 
    scatter3(ptsGood(:,1),ptsGood(:,2),ptsGood(:,3),'filled');
    scatter3(P(1),P(2),P(3),'blue','filled');
    hold off; view(3); xlabel('x'); ylabel ('y'); zlabel('z');
    title('Found Points');
subplot(1,2,2); scatter(u_raw,v_raw,'filled'); xlabel('u'); ylabel('v');
    title('Projected Points in UV Space');
%}
end