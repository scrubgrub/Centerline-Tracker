clear; clc; close all
% PARAMETERS
pts = 500;
sz1 = 224; sz2 = 224; sz3 = 224;

% Create ground truth XYZ
t = linspace(0,200,pts)';
offsetX = 112; offsetY = 112; offsetZ = 10;
x = offsetX+zeros(size(t)); y = offsetY+zeros(size(t)); z = offsetZ+t;
truthXYZ = [x y z];

sausage = zeros([sz1 sz2 sz3]);
coord_idx = sub2ind([sz1, sz2, sz3], round(x), round(y), round(z));
sausage(coord_idx) = 1;

radii = 4*[1 1.5 2 2.75 3 4 4.75 5 4.75 4 3 2.75 2 1.5 1 0.75]; % 9 different radii
for sli = 1:sz3
    rad = radii(mod(sli, length(radii)) + 1);
    se = strel('disk', rad);
    sausage(:,:,sli) = imdilate(squeeze(sausage(:,:,sli)), se);
end

volshow(sausage)

% Voxelize and identify unique coordinate points
%{
coord_idx = [round(x), round(y), round(z)];
coord_idx = unique(coord_idx, 'stable', 'rows');
coord_idx_lin = sub2ind([sz1, sz2, sz3], coord_idx(:,1), ...
    coord_idx(:,2), coord_idx(:,3));

theta = linspace(pi/6, 5*pi/6, pts);

% Create the sausage
sausage = extendedObjectMesh('sphere'); % Create first object out of loop
sausage = scale(sausage, sin(theta(1)));
sausage = translate(sausage, [coord_idx(1,2), ...
    coord_idx(1,1), coord_idx(1,3)]);
for pt = 2:size(coord_idx, 1) % Iteratively join future objects
    sphere = extendedObjectMesh('sphere');
    sphere = scale(sphere, sin(theta(pt)));
    sphere = translate(sphere, [coord_idx(pt,2), ...
        coord_idx(pt,1), coord_idx(pt,3)]);
    sausage = join(sausage, sphere);
end
show(sausage);

% Not functional - voxelization creates empty matrix
s = struct('vertices', sausage.Vertices, 'faces', sausage.Faces);
sausage = double(VOXELISE(1:sz1, 1:sz2, 1:sz3, s)); clear s
figure; hold on;
for sli = 1:sz3
    subplot(16, 28, sli);
    imshow(sausage(:,:,sli))
end

for sli = min(coord_idx(:,3)):max(coord_idx(:,3))
    sausage(:,:,sli) = imfill(squeeze(sausage(:,:,sli)));
end
%}

% Save Centerline
save("sausage.mat", "truthXYZ");
% Save Volume
niftiwrite(sausage, ...
    'C:\Users\jmath\OneDrive\Documents\MATLAB\MATLAB_SAVES\Sadleir\Tracking_Validation\simulateWires\sausage\new_sausage', ...
    "Compressed",true);