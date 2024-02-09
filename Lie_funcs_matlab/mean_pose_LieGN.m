% Optimization using Lie formulation. Gauss-Newton's method applied to the problem

function [T_op, cost_evol] = mean_pose_LieGN(T, weights, num_iters, verbosity)

num_poses = size(T,3);
if nargin < 4
    verbosity = false;
end
if nargin < 3
    num_iters = 5;
end
if nargin < 2
    weights = ones(num_poses,1);
end
weights = weights / sum(weights); % normalize weights to sum 1
assert(numel(weights) == num_poses);

% Initialize
% Use the middle pose as offset
idx_offset = ceil(num_poses/2);
T_op = T(:,:,idx_offset);
cost_evol = -ones(num_iters+1,1);
if verbosity
    cost_evol(1) = squared_distance_poses(T_op, T, weights);
    disp(['iter: ' num2str(0), '  cost = ' num2str(cost_evol(1))]);
end

% Iterate
for iter = 1:num_iters
    
    % Cast the problem into the tangent space (Lie Algebra)
    % Linearization. Analytical derivatives.
    epsilon = zeros(6,1);
    denom = 0;
    for k = 1:num_poses
        ww = weights(k)^2;
        epsilon = epsilon + ww * vee_se3( logm( T(:,:,k) / T_op ) );
        denom = denom + ww;
    end
    % Solve for the perturbation epsilon in the linear system of equations
    assert(denom > 0);
    epsilon = epsilon / denom;
    
    % Update "operating point"
    T_op = expm( hat_se3(epsilon) ) * T_op;
    
    if verbosity
        % Compute cost, to show how it evolves
        cost_evol(1+iter) = squared_distance_poses(T_op,T, weights);
        disp(['iter: ' num2str(iter), '  cost = ' num2str(cost_evol(1+iter)), '  norm(epsilon) = ' num2str(norm(epsilon))]);
    end
end
