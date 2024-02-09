% Point cloud registration using poses, with and without Lie Theory
%
% State Estimation for Robotics
% TU Berlin
% Guillermo Gallego

clear, clc, close ALL

disp('Point cloud registration in SE(3)');

num_points = 12;
x_true = randn(3,num_points); % 3D points
xh_true = [x_true; ones(1,num_points)]; % 3D points in homog coords

% Generate the ground truth pose
phi_exact = pi * randn(3,1);
rho_exact = 2 * randn(3,1);
xi_exact = [rho_exact; phi_exact];
T_exact = expm(hat_se3(xi_exact))

yh_true = T_exact * xh_true; % transformed (rot and trans) 3D points

%% Add noise
sigma = 0.1; % standard deviation of noise
noise = sigma * randn(4,num_points);
noise(4,:) = 0;

yh = yh_true + noise;

% Before alignment / registration
figure,
plot3(xh_true(1,:),xh_true(2,:),xh_true(3,:),'-ob')
hold on
plot3(yh(1,:),yh(2,:),yh(3,:),'-or')
axis vis3d
daspect([1 1 1])
grid
xlabel('X'), ylabel('Y'), zlabel('Z')
title('Before optimization')


%% Optimization in R^6 parameter space: search for the pose vector
% Yes, for this registration problem with poses there is a closed-form
% solution. Ignore it.
disp('-------- FMINUNC');
disp('Optimization in R^6 parameter space: search for the pose vector');

% Objective function: sum of square distances between points
costf = @(xi) cost_func_se3(xi,xh_true,yh);

% Optimize
options = optimoptions('fminunc','Display','iter', 'PlotFcn','optimplotfval');
r0 = zeros(6,1) % initial parameter vector in the optimization
%r0 = rotvec_exact + randn(6,1)
[r,fval,exitflag,output] = fminunc(costf,r0,options);
disp('best parameter vector:');
r

% Get pose matrix from best pose vector
T_param_space_search = expm(hat_se3(r));
disp('Distance between exact pose and estimated pose');
delta = norm(vee_se3(logm(T_exact / T_param_space_search)))


%% 
disp('-------- LSQNONLIN');
disp('Optimization in R^6 parameter space: search for the pose vector');

% Objective function: sum of square distances between points
res_func = @(xi) residual_func_se3(xi,xh_true,yh);

% Optimize
options = optimoptions('lsqnonlin','Display','iter', 'PlotFcn','optimplotresnorm');
%options.Algorithm = 'trust-region-reflective';
options.Algorithm = 'levenberg-marquardt';
[r,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(res_func,r0,[],[],options);
disp('best pose vector (lsqnonlin):');
r
h_optim_plot = gcf;

% Get pose matrix from best pose vector
T_param_space_search = expm(hat_se3(r));
disp('Distance between exact pose and estimated pose');
delta = norm(vee_se3(logm(T_exact / T_param_space_search)))

% Covariance
if sigma > 0
    covar_mat = inv( (1/(sigma^2))*(jacobian'*jacobian) );
    sqrt(diag(covar_mat))
end

%% Optimization using Lie formulation. Newton's method adapted to poses
disp('--------');
disp('Optimization using Lie-sensitive perturbations');
num_iters = 10;
cost = -ones(num_iters+1,1);
T_op = expm(hat_se3(r0));
cost(1) = squared_distance(T_op*xh_true,yh);
disp(['iter: ' num2str(0), '  cost = ' num2str(cost(1))]);

% Iterate
for iter = 1:num_iters

    % Compute the linear system Hessian * perturbation = - gradient 
    A = zeros(6,6);
    b = zeros(6,1);
    for j=1:num_points
        zhj = T_op * xh_true(:,j);
        Zj = circledot(zhj);
        A = A + Zj'*Zj;
        b = b + Zj'*(yh(:,j) - zhj);
    end
    % Solve for the perturbation epsilon in the linear system of equations
    epsilon = A \ b;
    
    % Update "operating point", pose matrix T_op
    T_op = expm(hat_se3(epsilon)) * T_op;
    
    % Compute cost, to show how it evolves
    cost(1+iter) = squared_distance(T_op*xh_true,yh); % compute new cost
    disp(['iter: ' num2str(iter), '  cost = ' num2str(cost(1+iter)), '  norm(epsilon) = ' num2str(norm(epsilon))]);
    
    % Display evolution of the optimization process
    yh_pred = T_op*xh_true;
    figure(100),
    plot3(yh(1,:),yh(2,:),yh(3,:),'-or')
    hold on
    plot3(yh_pred(1,:),yh_pred(2,:),yh_pred(3,:),'-ob')
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

figure(h_optim_plot)
hold on
plot(0:num_iters, cost ,'.-')
grid on
xlabel('iteration')
title('Lie-sensitive perturbations: Evolution of cost value')
if sigma < 1e-5
    h = gca; h.YScale = 'log'
end

disp('Distance between exact and estimated pose');
delta_param_Lie = norm(vee_se3(logm(T_exact / T_op)))

% The solutions obtained by the two methods are numerically very close:
disp('Distance between both estimated poses');
delta = norm(vee_se3(logm(T_param_space_search / T_op)))
