% Examples of how to generate random rotations using Lie Theory
%
% State Estimation for Robotics
% TU Berlin
% Guillermo Gallego

clear, clc, close ALL

% A 3D point (to be rotated)
x = [1,0,0]';

% Mean rotation
%C_mean = expm(pi*hat_so3(randn(3,1)));
C_mean = eye(3);

% Random samples in the tangent space (parameter space or Lie Algebra)
sigma = .1; % standard deviation of noise
num_samples = 500;
e_samples = randn(3,num_samples);
%e_samples = sigma * e_samples; % isotropic
e_samples = diag([sigma,sigma,2*sigma])* e_samples; % anisotropic

% Comute the random rotation and apply it to the point
y = zeros(3, num_samples);
for k = 1:num_samples
    y(:,k) = expm(hat_so3(e_samples(:,k))) * C_mean * x;
end

ym = C_mean * x; % point rotated using mean rotation

% Plot rotated point(s)
figure,
sphere
axis equal
hold on, 
plot3(ym(1),ym(2),ym(3),'*k')
plot3(y(1,:),y(2,:),y(3,:),'.r')
axis vis3d
daspect([1 1 1])
grid on
xlabel('X'), ylabel('Y'), zlabel('Z')
