% Smoothing of poses using Lie Theory over a sliding window
%
% State Estimation for Robotics
% TU Berlin
% Guillermo Gallego

clc, clear, close ALL

addpath('../Lie_funcs_matlab');

% profile on

%% Loading data
num_poses = 2000;
idx_first_pose = 6000;
idx_last_pose = idx_first_pose + num_poses;

% Real data
%groundtruth = importfile('data/shapes_rotation/groundtruth.txt', idx_first_pose, idx_last_pose);
groundtruth = importfile('data/poster/groundtruth.txt', idx_first_pose, idx_last_pose);
q = groundtruth(:,[8,5:7]);
RotMats = quat2rotm(q);
tvecs = groundtruth(:,[2,3,4])';
PoseMats = nan*ones(4,4,num_poses);
for k = 1:num_poses
    PoseMats(:,:,k) = [RotMats(:,:,k), tvecs(:,k); 0 0 0 1];
end


%% Filtering the data
N = 11;
N_half = ceil(N/2); % N odd
idx_valid = N_half:(num_poses-N_half);

% If plain average
PoseMats_smooth = smooth_poses(PoseMats, N);
posevecs_smooth_rect = PoseMats2PoseVecs(PoseMats_smooth(:,:,idx_valid));

% If weighted (e.g., Gaussian) average. See window(@WNAME,N) function
weights = gausswin(N);
% weights = hanning(N);
PoseMats_smooth = smooth_poses(PoseMats, N, weights);
posevecs_smooth_gauss = PoseMats2PoseVecs(PoseMats_smooth(:,:,idx_valid));

posevecs_orig = PoseMats2PoseVecs(PoseMats);


% profile viewer
% The logm function is consuming 80% of the processing. -> We should use a
% faster, specific one for this Lie Group.

%% Visualize DOFs
figure,
plot(1:size(posevecs_orig,2), posevecs_orig(1:3,:)')
hold on,

plot(N_half+(0:size(posevecs_smooth_rect,2) -1),posevecs_smooth_rect(1:3,:)'),
plot(N_half+(0:size(posevecs_smooth_gauss,2)-1),posevecs_smooth_gauss(1:3,:)','--')
title('exp coords - rho')
legend('X','Y','Z', 'X smooth rect','Y smooth rect','Z smooth rect', ...
    'X smooth Gauss','Y smooth Gauss','Z smooth Gauss')

figure,
plot(1:size(posevecs_orig,2), posevecs_orig(4:6,:)')
hold on,
plot(N_half+(0:size(posevecs_smooth_rect,2) -1),posevecs_smooth_rect(4:6,:)'), 
plot(N_half+(0:size(posevecs_smooth_gauss,2)-1),posevecs_smooth_gauss(4:6,:)','--')
title('exp coords - phi (rotation)')
legend('X','Y','Z', 'X smooth rect','Y smooth rect','Z smooth rect', ...
    'X smooth Gauss','Y smooth Gauss','Z smooth Gauss')
