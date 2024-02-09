% Returns a real number so that it can be optimized using fminunc

function res = cost_func_se3(v, xh,yh)
T = expm(hat_se3(v)); % pose homog matrix
yh_predicted = T*xh;
res = squared_distance(yh(1:3,:), yh_predicted(1:3,:));
