Author: Warren Boschen & Jade Lariviere, February 20 2026

Sets:
    ASU: 2, 4, 22, 23, 26, 47, 51, 54, 64, 70, 73, 75, 76, 86, 92, 97, 100, 1001
    UF: 1, 2, 4, 5, 6, 7
    Test: branching (synthetic and real), loop, spiral, sausage, wavy

Each subject folder for each set (ASU, UF, Test) contains the same files:
    SUB###_WIRE.mat - Centerline produced by the WireTrack script (ORTHO).
    SUB###_WIRE_VMTK.mat - Centerline produced by VMTK in 3D Slicer.
    *.json - Direct output of centerline produced by VMTK. Read into MATLAB using read_vmtk.m and produces SUB###_WIRE_VMTK.mat.
    *.nii.gz - The binary masks of the listed wire (F3, F4, OZ for all except ASU1001, who has FPZ, OZ, T7, T8). File names are structured differently between sets.