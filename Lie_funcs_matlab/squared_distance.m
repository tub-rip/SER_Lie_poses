function res = squared_distance(x,y)
num_points = size(x,2);
res = sum(sum((x - y).^2)) / num_points;
