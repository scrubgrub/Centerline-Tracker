% Author: Warren Boschen
% Date: June 26, 2025
%   Edited by Jade Lariviere (Feb. 20, 2026)
%{
Perform quantitative analysis of centerline tracking compared to the
ground truth.
===========================================================================
USAGE: Hit Run at the top of the screen.

INPUTS: Change PARAMETERS section if necessary to configure to data, and 
lines 60-65 ("LOAD DATA HERE!!") then Run. UF, ASU and test data had slight
differences between labeling originally, hence the hardcoded lines specific
to each dataset.

OUTPUT: In the Command Window, the following results are listed:
    - Mean & standard deviation of absolute difference between ground truth
        and tracked centerline.
    - Mean & standard deviation of percent difference between ground truth
        and tracked centerline.
    - RCAAEF Overlap between centerlines. Only works correctly if has TRUTH.
    - Sensitivity, Precision, Accuracy, Overlap (Dice coefficient).
    - Ratio of maximally-inscribed sphere reconstructed vol to original vol.

DEPENDENCIES: this code requires the VOXELISE function and its subfunctions
created by Adam A. on MATLAB File Exchange.
    (www.mathworks.com/matlabcentral/fileexchange/27390-mesh-voxelisation)
===========================================================================
%}

%% PARAMETERS
clear; clc; close all

subj= [2 4 22 26 47 51 54 64 70 73 75 76 86 92 97 100 1001];
    geometry = ["F3" "F4" "OZ"];
%subj = [1001]; geometry = ["FPZ" "OZ" "T7" "T8"];
% subj = "test"; geometry = ["wavy" "spiral" "loop" "sausage"];

set_name = "ASU"; 

hasTruth = false; % does the volume have a ground truth centerline?
shouldPrint = true; % print and plot analysis results?
%% Calculating Metrics
filePath = matlab.desktop.editor.getActiveFilename;
    [filePath,~,~] = fileparts(filePath);
    addpath(filePath);
tableList = cell([length(subj) , length(geometry)]);

ID = sprintfc("%03u",subj);
    % ID = sprintfc("%s",subj); % for test shapes

I = length(geometry); % preallocate due to parfor loop
ptDiffs = cell([length(geometry) 2]);
tic;
for n = 1:length(subj)
    Filepath = "C:\Users\###\SUB"+ID{n};
    addpath(genpath(Filepath));
    for m = 1:I
        % LOAD DATA HERE!!
        if hasTruth, TRUTH = load("true_"+geometry(m)+".mat").truthXYZ; end %#ok<UNRCH>
        % -------------------------------------------------------------------------
        wire_mask = niftiread("FULL_MODEL_"+ID{n}+"_MREIT_"+geometry(m)+" Wire.nii.gz"); % ASU
            wire_mask = wire_mask/max(wire_mask,[],"all"); % for ASU
        % wire_mask = niftiread("wire_"+geometry(m)+"_closed_reg.nii.gz"); % for UF
        % wire_mask = niftiread(geometry(m)); % for test shapes
        % -------------------------------------------------------------------------
        % remove extraneous CC in image, keep largest one
        wire_CCs = labelmatrix(bwconncomp(wire_mask));
        biggest_CC = mode(wire_CCs(wire_CCs ~= 0),"all");
        wire_mask(wire_CCs ~= biggest_CC) = 0;

        % load centerlines
        VMTK = load(sprintf("SUB%s_%s_VMTK.mat",ID{n},geometry(m))).vmtkCentroids;
        ORTHO = load(sprintf("newSUB%s_%s.mat",ID{n},geometry(m))).centerlines.Coordinates;
            % VMTK = load(geometry(m)+"_vmtk.mat").vmtkCentroids; % test shape
            % ORTHO = load(geometry(m)+"_ortho.mat").centerlines.Coordinates; % test shape
        % Comparisons =====================================================
        % VMTK ------------------------------------------------------------
        fprintf('\n========== VMTK RESULTS ==========\n%s_%s\n',ID{n},geometry(m));
        if hasTruth
           [VMTKdiffs,~] = calcDistances(VMTK, TRUTH, shouldPrint); %#ok<UNRCH>
           VMTKoverlap = Overlap(VMTK, TRUTH, wire_mask, shouldPrint);
        else
            [VMTKdiffs,~] = calcDistances(VMTK, nan, shouldPrint);
            VMTKoverlap = Overlap(VMTK, nan, wire_mask, shouldPrint);
        end 
        [~,VMTKstats] = maxSpheres(VMTK, wire_mask, shouldPrint,'blue',set_name+ID{n}+geometry(m));
        VMTKpts = size(VMTK,1); fprintf('VMTK numPts = %d\n',VMTKpts);
        
        % Ortho -----------------------------------------------------------
        fprintf('\n========== ORTHO RESULTS ==========\n%s_%s\n',ID{n},geometry(m));
        if hasTruth
            [WIREdiffs,~] = calcDistances(ORTHO, TRUTH, shouldPrint); %#ok<UNRCH>
            WIREoverlap = Overlap(ORTHO, TRUTH, wire_mask, shouldPrint);
        else
            [WIREdiffs,~] = calcDistances(ORTHO, nan, shouldPrint); 
            WIREoverlap = Overlap(ORTHO, nan, wire_mask, shouldPrint);
        end
        [~,WIREstats] = maxSpheres(ORTHO, wire_mask, shouldPrint, 'red', set_name+ID{n}+geometry(m));
        ORTHOpts = size(ORTHO,1); fprintf('ORTHO numPts = %d\n',ORTHOpts);
        
        % make table for export -------------------------------------------
        % VMTK input first, followed by ORTHO method
        tableList{n,m} = makeTable({VMTKdiffs; WIREdiffs},{VMTKstats; WIREstats}, ...
            {VMTKoverlap; WIREoverlap}, {VMTKpts; ORTHOpts});

        % run stats analysis (test shape only) ----------------------------
        % [h,p,ci,stats] = runTtest(TRUTH, ORTHO, VMTK);
        % fprintf("ORTHO vs VMTK (%s): h = %d, p = %d\n",geometry(i),h,p)
    end
end
toc;

%% Export Stat Tables
% exported as one file with sheets for each geometry
for n = 1:size(tableList, 1)
    for m = 1:length(geometry)
        writetable(tableList{n,m},strcat("statTable_",set_name,".xlsx"), ...
            "Sheet",ID{n},"WriteRowNames",true,"Range", ...
            sprintf("A%u:P%u",(4*(m-1)+1),((4*(m-1)+3))));
    end
end

%% Functions
function [stats, diff_norm] = calcDistances(Est, True, doPrint)
    [diff, idx] = distsToTruth(Est, True);
    stats.std_diff = std(diff); stats.mean_diff = mean(diff);
    pdiff = 2*(abs(Est - True(idx,:))./(True(idx,:)))*100;
    stats.mean_pdiff = mean(mean(pdiff,1)); stats.std_pdiff = std(pdiff, 0, "all");
    [stats.max_diff, max_idx] = max(diff);
    stats.max_pdiff = sqrt(pdiff(max_idx,1).^2 + pdiff(max_idx,2).^2 ...
                + pdiff(max_idx,3).^2);
    diff_norm = diff./size(Est,1); % normalized to diffs per point
    stats.mean_diff_norm = mean(diff_norm);
    stats.std_diff_norm = std(diff_norm);
    if doPrint
        fprintf('Mean distance to ground truth: %0.4f ± %0.4f (SD).\n', ...
            stats.mean_diff, stats.std_diff);
        fprintf('Mean percent difference is %0.4f%% ± %0.4f%% (SD). \n', ...
            stats.mean_pdiff, stats.std_pdiff);
        fprintf('Max distance to ground truth is %0.4f. \n', stats.max_diff);
        fprintf('Max percent difference is %0.4f%%. \n', stats.max_pdiff);
    end
end

function [h, p, ci, stats] = runTtest(True, ORTHO, VMTK)
    [ORTHOdiff, ~] = distsToTruth(ORTHO, True);
    ORTHOdiff_norm = ORTHOdiff./size(ORTHO,1); % normalized to diffs per point
    [VMTKdiff, ~] = distsToTruth(VMTK, True);
    VMTKdiff_norm = VMTKdiff./size(VMTK,1); % normalized to diffs per point
    [h, p, ci, stats] = ttest2(ORTHOdiff, VMTKdiff, ...
        "Alpha",0.05,"Tail","both");
end

function [maxRad,stats] = maxSpheres(Est, Vol, plotFig, plotColor, Name)
    [sz1, sz2, sz3] = size(Vol);
    dist_transform = double(bwdist(~Vol));
        dist_transform(dist_transform==0) = 1e-3; % ensure no zeros
    
    % Voxelize coordinates from Bezier curve
    coord_idx = [round(Est(:,1)), round(Est(:,2)), round(Est(:,3))];
    coord_idx = unique(coord_idx, 'stable', 'rows');
    coord_idx_lin = sub2ind([sz1, sz2, sz3], coord_idx(:,1), ...
                    coord_idx(:,2), coord_idx(:,3));
    
    % Create extendedObjectMesh of maximally inscribed spheres
    tube = extendedObjectMesh('sphere'); % Create first object out of loop
    tube = scale(tube, dist_transform(coord_idx_lin(1)));
    tube = translate(tube, [coord_idx(1,2), ...
                     coord_idx(1,1), coord_idx(1,3)]);
    maxRad = dist_transform(coord_idx_lin(1));
    for pt = 2:size(coord_idx, 1) % Iteratively join future objects
        if dist_transform(coord_idx_lin(pt)) > 0 % Double check
            sphere = extendedObjectMesh('sphere');
            sphere = scale(sphere, dist_transform(coord_idx_lin(pt)));
            sphere = translate(sphere, [coord_idx(pt,2), ...
                     coord_idx(pt,1), coord_idx(pt,3)]);
            tube = join(tube, sphere);
            
            % Search for largest radius of maximally inscribed sphere
            if dist_transform(coord_idx_lin(pt) > maxRad)
                maxRad = dist_transform(coord_idx_lin(pt));
            end
        end
    end

    % Calculate volume (better way to do this?)
    [F, V] = isosurface(Vol);
    [~, volTrue] = convhulln(V);
    [~, volExperimental] = convhulln(tube.Vertices);
    stats.volRatio = volExperimental/volTrue;
    
    % Voxelize complete tube structure, multiple passes
    t = struct('vertices', tube.Vertices, 'faces', tube.Faces);
    grid1 = VOXELISE(1:sz1, 1:sz2, 1:sz3, t, 'xy');
    grid2 = VOXELISE(1:sz1, 1:sz2, 1:sz3, t, 'yz');
    grid3 = VOXELISE(1:sz1, 1:sz2, 1:sz3, t, 'zx');
    tubeGrid = (grid1 + grid2 + grid3) > 1;
    
    % Perform morphological operation cleanup
    SE = strel('sphere', 2);
    tubeGrid = imclose(tubeGrid,SE);

    % Ensure that ifsmember does not include values outside of true positive
    % DOI: 10.1109/TIP.2018.2862346
    intersection = sum(tubeGrid.*permute(Vol,[2 1 3]), "all");
    stats.sensitivity = intersection/sum(Vol, "all")*100; 
    stats.precision = intersection/sum(tubeGrid, "all")*100; 
    stats.accuracy = ((stats.sensitivity + stats.precision)/2);
    % Overlap here is equivalent to Dice coefficient
    stats.overlap = ((2*intersection)/(sum(tubeGrid, "all") + ...
        sum(Vol, "all")))*100;
    
    % Plot
    if plotFig
        fprintf('Sensitivity: %0.2f%%\n', stats.sensitivity)
        fprintf('Precision: %0.2f%%\n', stats.precision)
        fprintf('Accuracy: %0.2f%%\n', stats.accuracy)
        fprintf('Overlap: %0.2f%%\n', stats.overlap)

        figure; hold on;
        fprintf('Ratio of tracked to true volume: %0.4f \n', stats.volRatio);
        patch('Faces', F, 'Vertices', V, 'FaceColor', 'green', ...
              'FaceAlpha', 0.3, 'EdgeColor', 'none');
        patch('Faces', tube.Faces, 'Vertices', tube.Vertices, 'FaceColor', ...
              plotColor, 'FaceAlpha', 0.7, 'EdgeColor', 'none');
        view(100,-10); axis equal; hold off; axis off; axis ij;
        subtitle(Name);
    end
end

function overlap = Overlap(Est, True, Vol, doPrint)
    % https://doi.org/10.3390/jimaging9120268
    % Perfect overlap is defined as the case when "all extracted centerline
    % points are located within one vessel radius of all reference
    % centerline points, and having no extraneous points."
     
    TPM = 0; % True positives of centerline
    TPR = size(True, 1); % True positives of reference
    FN = 0; % False negative
    FP = 0; % False positives
    % NOTE: A point of the reference standard is marked as TPR if
    % the distance to at least one of the connected points on the evaluated
    % centerline is less than the annotated radius and FN otherwise.
    % Since the volume was created FROM the known centerline, all points on
    % the "True" line are TPR and none are FN.
    
    radius = maxSpheres(Est, Vol, false, "");
    [minDist, ~] = distsToTruth(Est, True);
    
    for pt = 1:size(Est, 1)
        if minDist(pt) < radius
            TPM = TPM + 1;
        else
            FP = FP+1;
        end
    end

    overlap = (TPM + TPR)/(TPM + TPR + FN + FP);
    if doPrint
        fprintf('RCAAEF Overlap between centerlines: %0.4f\n',overlap)
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

function statTable = makeTable(diffs,stats,overlap,pts)
    % helper function to compile stats from both for export
    diffs = [diffs{1}; diffs{2}]; stats = [stats{1}; stats{2}]; 
    overlap = [overlap{1}; overlap{2}]; pts = [pts{1}; pts{2}];

    diffsTable = struct2table(diffs);
    statsTable = struct2table(stats);
    overlapTable = table(overlap,'VariableNames',"RCAAEF overlap");
    pointsTable = table(pts, 'VariableNames',"numPoints");
    statTable = [diffsTable statsTable overlapTable pointsTable];
    statTable.Properties.VariableUnits = ...
            {'voxels','voxels','%','%','voxels','%','voxels/pt','voxels/pt','','%','%','%','%','',''};
    statTable.Properties.RowNames = {'VMTK','ORTHO'};
end