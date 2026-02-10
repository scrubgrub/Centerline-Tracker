%% create test 3D volume
% Jade Lariviere || last modified Mar. 26, 2025
clc; clear; close all;
% SET UP ==================================================================
% name = inputdlg('Shape Name','Input shape name.'); % get wire name
% create parameterized centerline coordinates for volume
% PARAMETERS FOR FUNCTION
t = linspace(-5,5,1500); % how much of truth centerline to get
% spiral = 0,25,1500
% wavey = 0,10,1500
% loop = -5,5,1500
% spiral: [112,112,10]; wavey: [112,112,10]; loop: [112,112,100];
offsetX = 112; offsetY = 112; offsetZ = 10;
% define function to use here: makeWavy, makeSpiral, makeLoop
[x,y,z] = makeLoop(t);
% RESULTS =================================================================
% center shape in positive window
X = x + offsetX; Y = y + offsetY; Z = z + offsetZ;
truthXYZ = [X' Y' Z'];
% convert into a 3D volume matrix -----------------------------------------
vol = zeros([224 224 224]);
coord_idx = sub2ind(size(vol),round(X),round(Y),round(Z));
vol(coord_idx) = 1;
% dilate to produce tube-like volume
se = strel('sphere',4);
vol = imdilate(vol,se);
% plot result to double-check
figure(1); % only shows mathematical center
plot3(X,Y,Z); xlabel('x'); ylabel('y'); zlabel('z');
volshow(vol); % test volume

% save centerline
save("sausage.mat", "truthXYZ");
% save volume
vol_nii = make_nii(vol);
save_nii(vol_nii, 'sausage.nii');
%% FUNCTIONS LIST =========================================================
% helper functions --------------------------------------------------------
function c_list = loadCentroids(path,name)
% helper function to load centroids from file
c = load(strcat(path,name));
c = struct2cell(c); % extract first entry
c_list = c{1};
if max(c_list,[],'all') < 1 % quick check if coords are < 1e-3
    c_list = c_list./1e-3; % scale back to normal
end
end
function [minDist, idxMin] = distsToTruth(EstXYZ,TruthXYZ)
% helper function to calculate closest dists to Truth + index
% create storage arrays for data
minDist = nan([length(EstXYZ) 1]);
idxMin = minDist; % equal size as minDist
for i = 1:length(EstXYZ) % for all points in estimate
    c = EstXYZ(i,:); c = repmat(c,[length(TruthXYZ),1]);
    [minDist(i),idxMin(i)] = min(vecnorm((TruthXYZ-c),2,2)); % calc dists
end
end
% shape functions ---------------------------------------------------------
function [x,y,z] = makeSpiral(t)
% PARAMETERS
r = 10; k = 2; % r = radius of spiral, k = z-axis stretch
x = r.*cos(t).*0.1.*t;
y = r.*sin(t).*0.1.*t;
z = k.*t;
end
function [x,y,z] = makeWavy(t)
% PARAMETERS
r = 10; k = 1; % r = radius, k = scaling
x = 2*r.*sin(t);
y = r.*(cos(t)+t);
z = t.^2./k;
end
function [x,y,z] = makeLoop(t)
% PARAMETERS
r = 15; k = 2; % r = radius, k = scaling
x = k*(t).^2;
y = r*cos(t)+(t*0.6);
z = 2*r*sin(t);
end