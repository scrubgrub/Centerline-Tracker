%% helper function to convert from UV space to XYZ space
% Jade Lariviere | last modified Aug. 4, 2025

function [XYZ] = UVtoXYZ(UV,Axes,PlaneOffset,ImageOffset)
% this function outputs an Nx3 array from Nx2 input UV array, using the 
% Nx3 PlaneOffset (x,y,z) and 1x2 ImageOffset (u,v) variables. PlaneOffset
% determines whether XYZ is an array of points or vectors. ================

% validate UV and PlaneOffset size compatibility
if all([(size(PlaneOffset,1) ~= 1), (size(PlaneOffset,1) ~= size(UV,1))], "all")
    error('UVtoXYZ(): PlaneOffset 1st dimension size must be 1 or equal to UV.');
end

% convert points
XYZ = (UV(:,1)-ImageOffset(1))*Axes.U + (UV(:,2)-ImageOffset(2))*Axes.V ...
    + PlaneOffset;

end