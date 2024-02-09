% Smoothing of rotations using Lie Theory over a sliding window
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
groundtruth = importfile('data/shapes_rotation/groundtruth.txt', idx_first_pose, idx_last_pose);
q = groundtruth(:,[8,5:7]);
RotMats = quat2rotm(q);

% % Random data
% RotMats = nan*ones(3,3,num_poses);
% RotMats(:,:,1) = eye(3);
% for k = 2:num_poses
%     RotMats(:,:,k) = RotMats(:,:,k-1) * expm(hat_so3(0.2*randn(3,1)));
%     %RotMats(:,:,k) = expm(hat_so3(randn(3,1)));
%     %norm(RotMats(:,:,k)'*RotMats(:,:,k) - eye(3),'fro')
% end


%% Filtering the data
N = 11;
N_half = ceil(N/2); % N odd
idx_valid = N_half:(num_poses-N_half);

% If plain average
RotMats_smooth = smooth_rotations(RotMats, N);
rotvecs_smooth_rect = RotMats2RotVecs(RotMats_smooth(:,:,idx_valid));

% If weighted (e.g., Gaussian) average. See window(@WNAME,N) function
weights = gausswin(N);
RotMats_smooth = smooth_rotations(RotMats, N, weights);
rotvecs_smooth_gauss = RotMats2RotVecs(RotMats_smooth(:,:,idx_valid));

rotvecs_orig = RotMats2RotVecs(RotMats);


% profile viewer
% The logm function is consuming 80% of the processing. -> We should use a
% faster, specific one for this Lie Group.

%% Visualize DOFs
figure,
plot(1:size(rotvecs_orig,2), rotvecs_orig')
% set(gca,'ColorOrderIndex',1)
hold on,

plot(N_half+(0:size(rotvecs_smooth_rect ,2)-1),rotvecs_smooth_rect')
plot(N_half+(0:size(rotvecs_smooth_gauss,2)-1),rotvecs_smooth_gauss','--')
title('Rotation DOFs')
legend('X','Y','Z', 'X smooth rect','Y smooth rect','Z smooth rect', ...
    'X smooth Gauss','Y smooth Gauss','Z smooth Gauss')
