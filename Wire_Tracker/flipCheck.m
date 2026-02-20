%% enforce continuous direction, prevent flipping on path
% Jade Lariviere | last modified Dec. 9, 2025

function [hasFlipped, P, N] = flipCheck(LastPoints)
% this function determines whether the tracker has flipped direction while
% marching, and returns a logical true/false based on the evaluation 
% criterion between two determined vectors. The array LastPoints and 
% LastVectors must have at least 3 points (3x3). It can also output the XYZ
% point P and vector N the tracker should continue from (AKA reverts to 
% first point + dir if flip detected) =====================================

% recall A*B=0 -> perpendicular, A*B=-1 -> opposite, A*B=1 -> identical
flip_thresh =   -0.2; % dot product threshold

numMemory = size(LastPoints,1);
if numMemory < 3 % validate size of input argument (on 1st dim)
    error('flipCheck(): LastPoints must have a minimum memory of 3 points.');
end

% create normalized vectors approximating curve of previous points
if mod(numMemory,2) == 1 % odd number points
    mid_idx = median(1:numMemory); % middle index
    vec_prev = LastPoints(mid_idx,:)-LastPoints(1,:);
    vec_next = LastPoints(end,:)-LastPoints(mid_idx,:);
else % even number of points
    mid_idx = floor(median(1:numMemory));
    vec_prev = LastPoints(mid_idx,:)-LastPoints(1,:);
    vec_next = LastPoints(end,:)-LastPoints(mid_idx+1,:);
end
% normalize vectors
vec_prev = vec_prev/norm(vec_prev);
vec_next = vec_next/norm(vec_prev);

% check for flip via criterion ============================================
vec_dot = dot(vec_prev,vec_next);
if vec_dot < flip_thresh
    hasFlipped = true; % negative dot product = flipped
    % revert P and N (not reusing N to prevent infinite flipCheck loop)
    P = LastPoints(1,:);
    N = vec_prev;
else 
    hasFlipped = false;
    % keep P as-is, return apparent direction tracker is going as N
    P = LastPoints(end,:); N = vec_next;
end

% plots for debugging =====================================================
%{
figure(16);
    quiver3(0,0,0,vector_start(1),vector_start(2), ...
        vector_start(3),"cyan"); hold on;
    quiver3(0,0,0,vector_end(1),vector_end(2),vector_end(3),"red");
        hold off; xlim([-1 1]); ylim([-1 1]); zlim([-1 1]);
        xlabel('x'); ylabel('y'); zlabel('z');
        title("First vs. Last Direction", ...
            sprintf("A*B = %.3f",vector_dot)); view(3);
        legend('Initial Direction','Last Direction');
%}
end