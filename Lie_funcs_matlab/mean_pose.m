% Compute (weighted) mean of a collection of poses T

function T_mean = mean_pose(T, T_offset, weights)

num_poses = size(T,3);

if nargin < 3
    % All poses are equally weighted
    weights = ones(num_poses,1);
    
    if nargin < 2
        % Use the middle pose as offset
        idx_offset = ceil(num_poses/2);
        T_offset = T(:,:,idx_offset);
    end
end
weights = weights / sum(weights); % normalize weights to sum 1
assert(numel(weights) == num_poses);

% Compute incremental poses with respect to an offset one
deltaT_vecs = zeros(6,num_poses);
for k = 1:num_poses
    deltaT_vecs(:,k) = vee_se3( logm( T(:,:,k) / T_offset ) );
end

% Averaging step
deltaT_vec_mean = sum(repmat(weights',6,1) .* deltaT_vecs,2);

% Restore the offset to compute the corresponding pose
deltaT_mean = expm( hat_se3(deltaT_vec_mean) );
T_mean = deltaT_mean * T_offset;
