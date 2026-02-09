%% interpolate from series of points on 3D curve to destination point
% Jade Lariviere | last modified Sep. 24, 2025

function NewBranch = interpBranch(Branch, Main, Param)
% this function appends an interpolated curve produced from one 3D line
% ([Nx3] Branch) in reference to another Main 3D line, using information 
% from Branch and Main like average spacing and orientation, to try 
% producing a natural interpolation. ======================================

min_vec = Param.PruneLength; % smallest # vectors possible from a Branch
momentum_weighting = flip(cumsum(1:min_vec)./sum(cumsum(1:min_vec)));

% we assume first coordinate in Branch array is intended to merge to the 
% closest point from the main branch due to findBranch & goTrack behavior.
branch_start = Branch(1,:);
[~,idx_destination] = min(vecnorm(branch_start-Main,2,2),[],'all');
main_destination = Main(idx_destination,:);

% calculate input branch characteristics
branch_numPt = size(Branch,1);
branch_length = sum(vecnorm(Branch(1:end-1,:)-Branch(2:end,:),2,2));
branch_dir = sum((Branch(1:min_vec,:)-Branch(2:min_vec+1,:)).*momentum_weighting');

% calculate output interpolated curve statistics
interp_length = norm(main_destination-branch_start);
interp_step = branch_length/branch_numPt; % approx interp step = curve step

% parameterize between start and end; prep for loop
t = linspace(0,1,round(interp_length/interp_step));
interp = [branch_start; nan(length(t),3)]; % store interp curve

for i = 1:length(t) % linearly blend start vector & destination vector
    N = (1-t(i)).*branch_dir + t(i).*(main_destination-interp(i,:));
        N = N/norm(N);
    interp(i+1,:) = interp(i,:) + N*interp_step; % travel one step
end

NewBranch = [main_destination; flip(interp(2:end,:),1); Branch];
end