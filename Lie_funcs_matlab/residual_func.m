% Returns a residual vector so that it can be optimized using lsqnonlin

function res = residual_func(rotv,x,y)
C = expm(hat_so3(rotv)); % rotation matrix
y_err = y - C*x;
res = y_err(:); % residual vector

% Optional: scale by the number of points
num_points = size(x,2);
res = res / sqrt(num_points);
