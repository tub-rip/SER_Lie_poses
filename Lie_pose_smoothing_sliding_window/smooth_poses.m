% PoseMats is 4x4xnum_poses
% N is the number of poses in a window (length of the Finite Impulse
% Response -FIR- filter).

function PoseMats_smooth = smooth_poses(PoseMats, N, weights)

if nargin < 3
    weights = ones(N,1);
end
weights = weights / sum(weights); % normalize weights to sum 1
assert(numel(weights) == N);

% Book memory for the output
PoseMats_smooth = nan*PoseMats;

% Loop through the poses
num_poses = size(PoseMats,3);
N_half = ceil(N/2); % N odd

for k = 1:num_poses
    % Check indices first
    if (k < N_half) || ((num_poses-k) < N_half)
        continue;
    end
    
    % Get subset of poses
    Poses_subset = PoseMats(:,:,k+(1:N)-N_half);
    
    % Compute mean pose
    Pose_mean = mean_pose(Poses_subset, weights); % 1 iter from middle pose is good enough
    %Pose_mean = mean_pose_LieGN(Poses_subset, weights); % iterating

    % Assign the output
    PoseMats_smooth(:,:,k) = Pose_mean;
end

end


% Function that computes and approximate mean rotation. Non-iterative.
function T_mean = mean_pose(Poses_subset, weights)

N = size(Poses_subset,3);

% Use the middle pose as offset
N_half = ceil(N/2); % N odd
Pose_off = Poses_subset(:,:,N_half);
% Compute differences between poses in the FIR window and the offset
delta_paramvecs = zeros(6,N);
for ii = 1:N
    delta_paramvecs(:,ii) = vee_se3(logm(Poses_subset(:,:,ii) / Pose_off));
end
% Averaging step
deta_paramvec_mean = sum(repmat(weights',6,1) .* delta_paramvecs,2);
% Restore the offset to compute the corresponding pose
T_mean = expm(hat_se3(deta_paramvec_mean)) * Pose_off;

end