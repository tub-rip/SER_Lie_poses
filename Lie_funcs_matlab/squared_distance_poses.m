% Return the sum of squared pose distances between T_op and each pose in T

function cost = squared_distance_poses(T_op, T, weights)

num_poses = size(T,3);
if nargin < 3
    weights = ones(num_poses,1);
end

cost = 0;
for k = 1:num_poses
    err_vec = weights(k) * vee_se3( logm( T(:,:,k) / T_op ) );
    cost = cost +  err_vec'*err_vec;
end
