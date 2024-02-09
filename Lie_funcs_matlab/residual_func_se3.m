% Returns a residual vector so that it can be optimized using lsqnonlin

function res = residual_func_se3(v, xh,yh)
T = expm(hat_se3(v)); % pose homog matrix
yh_predicted = T*xh;
y_err = yh(1:3,:) - yh_predicted(1:3,:); % Euclidean coords
res = y_err(:); % residual vector

% Optional: scale by the number of points
num_points = size(xh,2);
res = res / sqrt(num_points);
