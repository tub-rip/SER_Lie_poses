% Examples of rotation or pose averaging using Lie Theory
%
% State Estimation for Robotics
% TU Berlin
% Guillermo Gallego

clear, clc, close ALL

%% Generate data
test1_interpolation_poses_LieTheory

disp('Pose averaging using Lie Theory')
num_poses

%% Rotation averaging. 1 pass or iteration

% Compute incremental rotations with respect to an offset one
% idx_offset = 1;
idx_offset = round(num_poses/2);
Rot_offset = Rots(:,:,idx_offset);
deltaR_vecs = zeros(3,num_poses);
for k = 1:num_poses
    deltaR = Rots(:,:,k) * Rot_offset';
    deltaR_vecs(:,k) = vee_so3( logm(deltaR) );
end
% Compute the mean of the incremental rotations
deltaR_vec_mean = mean(deltaR_vecs,2);
% Restore the offset to compute the corresponding rotation
deltaR_mean = expm(hat_so3(deltaR_vec_mean));
R_mean = deltaR_mean * Rot_offset;

% figure, plot3(deltaR_vecs(1,:),deltaR_vecs(2,:),deltaR_vecs(3,:),'.')

% Measure the distance between the "mean rotation" and each rotation
dist = zeros(1,num_poses);
for k = 1:num_poses
    deltaR = R_mean * Rots(:,:,k)';
    dist(k) = norm(vee_so3( logm(deltaR) ));
end
figure('Color','white') 
plot(dist*(180/pi),'linewidth',2), grid; % even or odd number of rotations?
xlabel('rotation index (k)')
ylabel('degrees')
title('Distance between the "mean rotation" and each rotation')

% Does the result depend on the offset rotation?
% With 1 iteration (and noise) => yes.


%% Pose averaging. 1 pass or iteration
% Same steps as for rotations, but now calling the corresponding hat and
% vee operators (for SE3).

T_offset = T(:,:,idx_offset);
T_mean = mean_pose(T, T_offset);

% Measure the distance between the "mean pose" and each pose
dist = zeros(1,num_poses);
for k = 1:num_poses
    deltaT = T_mean * inv(T(:,:,k));
    dist(k) = norm(logm(deltaT),'fro');
end
figure('Color','white') 
plot(dist,'linewidth',2), grid
xlabel('pose index (k)')
title('Distance between the "mean pose" and each pose')

 
%% Visualize how the coordinate axes transform using the interpolated poses
figure(1),
for k = 1:1
    p = T_mean * xyz_homog;
    hold on,
    plot3(p(1,1:2),p(2,1:2),p(3,1:2),'r','linewidth',2)  % plot x axis in red
    plot3(p(1,[1,3]),p(2,[1,3]),p(3,[1,3]),'g','linewidth',2) % plot y axis in green
    plot3(p(1,[1,4]),p(2,[1,4]),p(3,[1,4]),'b','linewidth',2) % plot z axis in blue
    axis vis3d, axis equal, grid on
    drawnow
end
hold off



%% Same thing, but with noisy poses

% Add noise to the poses
sigma = 0.1; % std of pose noise
T_noisy = T;
for k = 1:num_poses
    perturbation = sigma * randn(6,1);
    T_noisy(:,:,k) = expm( hat_se3( perturbation ) ) * T(:,:,k);
end

% Plot noisy poses
figure(11),
for k = 1:num_poses
    p = T_noisy(:,:,k) * xyz_homog;
    hold on,
    plot3(p(1,1:2),p(2,1:2),p(3,1:2),'r')  % plot x axis in red
    plot3(p(1,[1,3]),p(2,[1,3]),p(3,[1,3]),'g') % plot y axis in green
    plot3(p(1,[1,4]),p(2,[1,4]),p(3,[1,4]),'b') % plot z axis in blue
    axis vis3d, axis equal, grid on
    drawnow
end
hold off
xlabel('X'),ylabel('Y'),zlabel('Z')
title('Transformation of the coordinate axes')

% Compute mean pose (given an offset one)
T_offset = T_noisy(:,:,idx_offset);
T_mean = mean_pose(T_noisy, T_offset);

% Measure the distance between the "mean pose" and each noisy pose
dist = zeros(1,num_poses);
for k = 1:num_poses
    deltaT = T_mean * inv(T_noisy(:,:,k));
    dist(k) = norm(logm(deltaT),'fro');
end
figure('Color','white') 
plot(dist,'linewidth',2), grid
xlabel('pose index (k)')
title('Distance between the "mean pose" and each noisy pose')


%% Weighted average (akin Gaussian window in signal processing)

% Centered, bell-shaped weights
weights = hanning(num_poses);

% % Skewed weights
% weights = zeros(1,num_poses);
% half_kernel_size = fix(num_poses/2);
% weights(half_kernel_size + (1:half_kernel_size)) = hanning(half_kernel_size);

% weights = weights/sum(weights);

figure, plot(weights)
xlabel('pose index (k)')
title('Weights used for averaging')

T_weighted_mean = mean_pose(T_noisy, T_offset, weights);

% Measure the distance between the "weighted mean pose" and each rotation
dist = zeros(1,num_poses);
for k = 1:num_poses
    deltaT = T_weighted_mean * inv(T_noisy(:,:,k));
    dist(k) = norm(logm(deltaT),'fro');
end
figure('Color','white') 
plot(dist,'linewidth',2), grid
xlabel('pose index (k)')
title('Distance between the "weighted mean pose" and each noisy pose')

T_mean
T_weighted_mean
norm(vee_se3( logm( T_mean * inv(T_weighted_mean) )))


%% Optimization using Lie formulation. Gauss-Newton's method applied to the problem
disp('--------');
disp('Optimization-based point of view using Lie-sensitive perturbations');

num_iters = 6;
verbosity = true;
[T_op, cost_evol] = mean_pose_LieGN(T_noisy, weights, num_iters,verbosity);
% The iteration may not always lead to a decrease of the objective function

norm(vee_se3( logm( T_mean * inv(T_op) )))


%% See also Weiszfeld algorithm
% https://johnwlambert.github.io/lie-groups/
