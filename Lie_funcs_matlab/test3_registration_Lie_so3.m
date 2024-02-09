% Point cloud registration using rotations, with and without Lie Theory
%
% State Estimation for Robotics
% TU Berlin
% Guillermo Gallego

clear, clc, close ALL

disp('Point cloud registration in SO(3)');

%% Generate data (3D points)
num_points = 60;

bDumbbell = false;
if bDumbbell
    Np = 10;
    x0 = rand(3,Np); % 3D points
    x1 = x0 + [3,0,0]';
    Rot_y = expm(hat_so3([0,pi,0]'));
    x_true = [x1, Rot_y*x1];
    num_points = size(x_true,2);
else
    x_true = rand(3,num_points); % 3D points
end

% Plot the point cloud
% figure,
% plot3(x_true(1,:),x_true(2,:),x_true(3,:),'-ob')
% axis vis3d
% daspect([1 1 1])
% grid
% xlabel('X'), ylabel('Y'), zlabel('Z')

% Generate the ground truth rotation
rotvec_exact = randn(3,1);
% rotvec_exact = rotvec_exact * (pi / norm(rotvec_exact));
% rotvec_exact = rotvec_exact + 0.01*randn(3,1);

% rotvec_exact = [0,pi*0.95,0]' % Rot angle close to the singularity

angle = norm(rotvec_exact) * (180/pi);
disp(['Rotation angle (exact): ' num2str(angle) ' degrees']);

C_exact = expm(hat_so3(rotvec_exact)) % Rotation matrix
y_true = C_exact * x_true; % rotated (transformed) 3D points

% if true
%     % Mess up data association
%     idx = 1:Np;
%     y_true = [y_true(:,Np+idx),y_true(:,idx)];
% end

%% Add noise
sigma = 0.01; % standard deviation of noise
noise = sigma * randn(3,num_points);
y = y_true + noise;

% % Create outliers
% num_outliers = 10;
% indices = randperm(num_points);
% idx_outliers = indices(1:num_outliers);
% y(:,idx_outliers) = randn(3,num_outliers);

% Plot 3D points before alignment / registration
figure,
plot3(x_true(1,:),x_true(2,:),x_true(3,:),'-ob')
hold on
plot3(y(1,:),y(2,:),y(3,:),'-or')
axis vis3d
daspect([1 1 1])
grid
xlabel('X'), ylabel('Y'), zlabel('Z')
title('Before optimization')


%% Optimization in R^3 parameter space: search for the rotation vector
% Yes, for this registration problem with rotations there is a closed-form
% solution. Ignore it.
disp('-------- FMINUNC');
disp('Optimization in R^3 parameter space: search for the rotation vector');

% Objective function: sum of square distances between points
costf = @(rotv) cost_func(rotv,x_true,y);

% Optimize
options = optimoptions('fminunc','Display','iter', 'PlotFcn','optimplotfval');
r0 = [0,0,0]' % initial parameter vector in the optimization
%r0 = rotvec_exact + randn(3,1)
[r,fval,exitflag,output] = fminunc(costf,r0,options);
disp('best rotation vector (fminunc):');
r

% Get rotation matrix from best rotation vector
C_param_space_search = expm(hat_so3(r));
disp('Distance between exact and estimated rotation');
delta_param_R3 = norm(vee_so3(logm(C_exact'*C_param_space_search)))

%% 
disp('-------- LSQNONLIN');
disp('Optimization in R^3 parameter space: search for the rotation vector');

% Objective function: sum of square distances between points
res_func = @(rotv) residual_func(rotv,x_true,y);

% Optimize
options = optimoptions('lsqnonlin','Display','iter', 'PlotFcn','optimplotresnorm');
%options.Algorithm = 'trust-region-reflective';
options.Algorithm = 'levenberg-marquardt';
[r,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(res_func,r0,[],[],options);
disp('best rotation vector (lsqnonlin):');
r
h_optim_plot = gcf;

% Get rotation matrix from best rotation vector
C_param_space_search = expm(hat_so3(r));
disp('Distance between exact and estimated rotation');
delta_param_R3 = norm(vee_so3(logm(C_exact'*C_param_space_search)))

% Covariance
if sigma > 0
    covar_mat = inv( (1/(sigma^2))*(jacobian'*jacobian) );
    sqrt(diag(covar_mat)) * (180/pi)
end


%% Optimization using Lie formulation. Newton's method adapted to rotations
disp('--------');
disp('Optimization using Lie-sensitive perturbations');
num_iters = 10;
cost = -ones(num_iters+1,1);
C_op = expm(hat_so3(r0));
cost(1) = squared_distance(C_op*x_true,y);
disp(['iter: ' num2str(0), '  cost = ' num2str(cost(1))]);

% Iterate
for iter = 1:num_iters

    % Compute the linear system Hessian * perturbation = - gradient 
    A = zeros(3,3);
    b = zeros(3,1);
    for j=1:num_points
        zj = C_op * x_true(:,j);
        Zj = hat_so3(zj);
        A = A + Zj'*Zj;
        b = b + Zj'*(zj - y(:,j));
    end
    % Solve for the perturbation epsilon in the linear system of equations
    epsilon = A \ b;
    
    % Update "operating point", rotation matrix C_op
    C_op = expm(hat_so3(epsilon)) * C_op;
    
    % Compute cost, to show how it evolves
    cost(1+iter) = squared_distance(C_op*x_true,y); % compute new cost
    disp(['iter: ' num2str(iter), '  cost = ' num2str(cost(1+iter)), '  norm(epsilon) = ' num2str(norm(epsilon))]);
    
    % Display evolution of the optimization process
    y_pred = C_op*x_true;
    figure(100),
    plot3(y(1,:),y(2,:),y(3,:),'-or')
    hold on
    plot3(y_pred(1,:),y_pred(2,:),y_pred(3,:),'-ob')
    hold off
    title(['iter: ' num2str(iter)])
    axis vis3d
    daspect([1 1 1])
    grid on
    xlabel('X'), ylabel('Y'), zlabel('Z')
    drawnow
    pause(0.5);

end
title('After optimization')

%figure, 
figure(h_optim_plot)
hold on
plot(0:num_iters, cost ,'.-')
grid on
xlabel('iteration')
title('Lie-sensitive perturbations: Evolution of cost value')
if sigma < 1e-5
    h = gca; h.YScale = 'log'
end

disp('Distance between exact and estimated rotation');
delta_param_Lie = norm(vee_so3(logm(C_exact'*C_op)))

% The solutions obtained by the two methods are numerically very close:
disp('Distance between both estimated rotations');
delta_rots = norm(vee_so3(logm(C_param_space_search'*C_op)))


%% Closed-form solver
disp('--------');
disp('Optimization using closed-form formula');

% Data matrix
p = x_true;
p_centroid = mean(x_true,2);
y_centroid = mean(y,2);

%Method 1 (using a for loop)
% W = zeros(3,3);
% for j=1:num_points
%     W = W + (y(:,j)-y_centroid) * (p(:,j)-p_centroid)';
% end

%Method 2: vectorize, matrix multiplication.
yc = y - y_centroid; % centered values
pc = p - p_centroid; % centered values
W = yc*pc'; % 3x3 matrix

% Decompose data matrix and build estimated rotation matrix
[U,D,V] = svd(W);
C_closed_form = U*diag([1 1 det(U)*det(V)])*V';

disp('Distance between exact and estimated rotation');
delta_param_Lie = norm(vee_so3(logm(C_exact'*C_closed_form)))

% The solutions obtained by the two methods are numerically very close:
disp('Distance between both estimated rotations');
delta_rots = norm(vee_so3(logm(C_closed_form'*C_op)))
