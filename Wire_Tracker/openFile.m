%% function that handles importing of files
% Jade Lariviere | last modified Nov. 26, 2025

function [I,Wire] = openFile
% function extracts relevant data from the file imported into the script, 
% including name, directory, and extension, and it will also store 
% details relevant to the tracking algorithm, including volume dimensions 
% and the boundary in 3D, in structure I. structure Wire is also returned
% which is the volume itself, normalized + binarized. =====================

% read file depending on extension; *** change .raw details as needed! ***
[I.name,I.path] = uigetfile({'*.raw;*.nii.gz;*.nii','Structural Data'}, ...
    'Select file with wire segmentation data.');
    [~,~,I.ext] = fileparts(I.name); % get input file info

switch I.ext
    case '.raw' % for .raw files ------------------------------------------
    Wire.raw = double(multibandread(strcat(I.path,'\',I.name), ...
    [224 224 224],'*uint8',0,'bsq','ieee-be')); % fill file info by hand!!!
        Wire.raw = Wire.raw./max(Wire.raw,[],'all'); % normalize
    otherwise % for .nii & .nii.gz files ----------------------------------
        wire_in = niftiread(strcat(I.path,'\',I.name));
            %wire_in(wire_in~=1) = 0; % TEST OPERATION FOR LABELED FILES
        Wire.raw = wire_in./max(wire_in,[],'all'); % normalize
        % binarize input just in case
        Wire.raw = double(imbinarize(wire_in, "adaptive"));
end
% extract image dimensions / boundary
I.dims = size(Wire.raw); Wire.dims = I.dims; 
I.bounds = bwperim(ones(I.dims),6);
fprintf(2,'*** Loaded %s !\n',I.name);
end