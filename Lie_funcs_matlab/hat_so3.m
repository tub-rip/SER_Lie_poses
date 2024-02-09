% Auxiliary "hat" (lift) function, from the vector space to the Lie algebra
%
% check: hat_so3([1,2,3])

function M = hat_so3(v)
M = [0,-v(3),v(2); v(3),0,-v(1); -v(2),v(1),0];
