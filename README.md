# Centerline-Tracker
Code for a projection-based algorithm for tracking tube-like volumes.

Please read the full LICENSE file for details regarding distribution, commericalization, and 
derivative works. Thank you!

This MATLAB code was developed in r2024b and requires the Image Processing Toolbox™, as various
processing operations in the computer vision rely on its functions. This is not a compiled 
executable and requires the MATLAB Compiler Runtime in the Editor to use in its current state.
Open the .m file "WireTrack_Script" and select "Run" (F5).

NOTE: when reading in .raw files, the openFile.m script must be adjusted manually with .raw file
header information, like dimensions, dataType, endianness, etc., to be correctly imported.
