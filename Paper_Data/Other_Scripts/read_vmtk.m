%% Process centerline data from 3D Slicer VMTK

% Read in and pre-process centerline data from 3D Slicer VMTK

% Author: Warren Boschen
%           edited by Jade Lariviere, Oct. 24, 2025
% Date: June 19, 2025
% Usage: Change the name of the output file on Line 25.
%        In the Command Window, call:
%        json = 'NAME_OF_FILE';
%        centroids = read_vmtk(json);
% Input: json - A .json file exported from 3D Slicer.
% Output: A .mat file with the centroid information, with chosen name.
% References: https://slicer.readthedocs.io/en/latest/user_guide/about.html
%% FUNCTION read_vmtk

function vmtkCentroids = read_vmtk(json)
    % Read file
    [path,name,~] = fileparts(json);
    fid = fopen(json);
    raw = fread(fid, inf);
    str = char(raw');
    fclose(fid);
    % Encodes data and return as a structure
    val = jsondecode(str);
    vmtkCentroids = zeros(val.markups.lastUsedControlPointNumber, 3);
    for i=1:val.markups.lastUsedControlPointNumber
        vmtkCentroids(i,:) = val.markups.controlPoints(i).position;
    end
    
    % Check to make sure there are no negative values
    if find(vmtkCentroids(:,1) < 0)
        vmtkCentroids(:,1) = -1*vmtkCentroids(:,1);
    end
    if find(vmtkCentroids(:,2) < 0)
        vmtkCentroids(:,2) = -1*vmtkCentroids(:,2);
    end
    if find(vmtkCentroids(:,3) < 0)
        vmtkCentroids(:,3) = -1*vmtkCentroids(:,3);
    end
    % add [1 1 1]: going from 0 index (python) to 1 index (MATLAB)
    vmtkCentroids = vmtkCentroids + [1 1 1];

    s = strcat(path,"\",name);
    save(s, "vmtkCentroids");
    fprintf(2,'Saved as %s!\n',s);
end