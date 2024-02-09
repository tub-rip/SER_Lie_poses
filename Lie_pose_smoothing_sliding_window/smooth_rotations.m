% RotMats is 3x3xnum_poses
% N is the number of poses in a window (length of the Finite Impulse
% Response -FIR- filter).

function RotMats_smooth = smooth_rotations(RotMats, N, weights)

if nargin < 3
    weights = ones(N,1);
end
weights = weights / sum(weights); % normalize weights to sum 1
assert(numel(weights) == N);

% Book memory for the output
RotMats_smooth = nan*RotMats;

% Loop through the poses
num_poses = size(RotMats,3);
N_half = ceil(N/2); % N odd

for k = 1:num_poses
    % Check indices first
    if (k < N_half) || ((num_poses-k) < N_half)
        continue;
    end
    
    % Get subset of rotations
    Rots_subset = RotMats(:,:,k+(1:N)-N_half);
    
    % Compute mean rotation
    Rot_mean = mean_rotation(Rots_subset, weights);
    
    % Assign the output
    RotMats_smooth(:,:,k) = Rot_mean;
end

end


% Function that computes and approximate mean rotation. Non-iterative.
function Rot_mean = mean_rotation(Rots_subset, weights)

N = size(Rots_subset,3);

% Use the middle pose as offset
N_half = ceil(N/2); % N odd
R_off = Rots_subset(:,:,N_half);
% Compute differences between poses in the FIR window and the offset
delta_rvecs = zeros(3,N);
for ii = 1:N
    delta_rvecs(:,ii) = vee_so3(logm(Rots_subset(:,:,ii)*R_off'));
end
% Averaging step
deta_rvec_mean = sum(repmat(weights',3,1) .* delta_rvecs,2);
% Restore the offset to compute the corresponding pose
Rot_mean = expm(hat_so3(deta_rvec_mean)) * R_off;

end