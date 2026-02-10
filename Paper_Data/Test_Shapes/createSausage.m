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

% Save Centerline
save("sausage.mat", "truthXYZ");
% Save Volume
niftiwrite(sausage, ...
    'C:\Users\jmath\OneDrive\Documents\MATLAB\MATLAB_SAVES\Sadleir\Tracking_Validation\simulateWires\sausage\new_sausage', ...
    "Compressed",true);