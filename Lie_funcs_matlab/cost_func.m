% Returns a real number so that it can be optimized using fminunc

function res = cost_func(rotv,x,y)
C = expm(hat_so3(rotv)); % rotation matrix
res = squared_distance(y, C*x);
