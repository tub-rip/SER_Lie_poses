% vee_se3 Compute the pose vector corresponding to a Lie Algebra matrix for a pose

function xi = vee_se3(V)
xi = [V(1:3, 4); vee_so3( V(1:3, 1:3) )];
