% Examples of Interpolation of poses using Lie Theory
%
% State Estimation for Robotics
% TU Berlin
% Guillermo Gallego

clear, clc, close ALL

disp('Interpolation of poses using Lie Theory')

% Choose random exponential coordinates of poses (6x1 parameter vectors)
a = randn(6,1);
b = randn(6,1);
b(1:3) = b(1:3)*10;

a =[0.4889
    1.0347
    0.7269
   -0.3034
    0.2939
   -0.7873]

b =[8.8840
  -11.4707
  -10.6887
   -0.8095
   -2.9443
    1.4384]
    
% Start and End poses (obtained using the exponential map)
Ta = expm( hat_se3(a) )
Tb = expm( hat_se3(b) )


%% Compute interpolated poses using Lie Theory, in SE(3)
num_poses = 31;
u = linspace(0,1,num_poses);
T = zeros(4,4,num_poses);
for k = 1:num_poses
    T(:,:,k) = expm(u(k) * logm(Tb * inv(Ta))) * Ta;
    
    % % Other ways to obtain the interpolated poses:
    % T2 = Ta * expm(u(k) * logm(inv(Ta) * Tb));
    % norm(T2 - T(:,:,k),'fro') % is zero (machine precision)
    % T2 = Tb * expm((1-u(k)) * logm(inv(Tb) * Ta));
    % norm(T2 - T(:,:,k),'fro') % is zero (machine precision)
end
% T  % print poses


%% Visualize how the coordinate axes transform using the interpolated poses
xyz = [0 0 0; 1 0 0; 0 1 0; 0 0 1]';
xyz_homog = [xyz; ones(1,4)];

figure(1),
for k = 1:num_poses
    p = T(:,:,k)*xyz_homog;
    hold on,
    plot3(p(1,1:2),p(2,1:2),p(3,1:2),'r')  % plot x axis in red
    plot3(p(1,[1,3]),p(2,[1,3]),p(3,[1,3]),'g') % plot y axis in green
    plot3(p(1,[1,4]),p(2,[1,4]),p(3,[1,4]),'b') % plot z axis in blue
    axis vis3d, axis equal, grid on
    drawnow
end
hold off
xlabel('X'),ylabel('Y'),zlabel('Z')
title('Transformation of the coordinate axes')


%% Visualuze how the vertices of a cube transform using the interpolated poses
vert_xy = [1 1; 1 -1; -1 -1; -1, 1]';
xyz_cube = [vert_xy, vert_xy; ones(1,4), -ones(1,4)];
xyz_cube_homog = [xyz_cube; ones(1,8)];

% Plot the 8 vertices of the cube: the top square (z=1) and the bottom square (z=-1)
cmap = colormap;
idx = round(linspace(1,size(cmap,1),num_poses));
idx_cube = [1:4,1,5:8,5];
figure(2),
for k = 1:num_poses
    p = T(:,:,k)*xyz_cube_homog;
    hold on,
    plot3(p(1,idx_cube),p(2,idx_cube),p(3,idx_cube),'.-','Color',cmap(idx(k),:))
    axis vis3d, axis equal, grid on
    drawnow
end
hold off
xlabel('X'),ylabel('Y'),zlabel('Z')


%% Comparison against interpolation in the Cartesian product SO(3) x R^3 
Rot_a = Ta(1:3,1:3);
Rot_b = Tb(1:3,1:3);
transl_a = Ta(1:3,4);
transl_b = Tb(1:3,4);

% Compute interpolated rotations and translations, separately
Rots = zeros(3,3,num_poses);
transls = zeros(3,num_poses);
for k = 1:num_poses
    Rots(:,:,k) = expm(u(k) * logm(Rot_b * (Rot_a'))) * Rot_a;
    transls(:,k) = (1-u(k))*transl_a + u(k)*transl_b; % they span a 3D line
end

% Check the difference between the poses interpolated in SE(3) vs in SO(3)xR^3
% for k = 1:num_poses
%     % The rotational part is the same
%     % norm(Rots(:,:,k) - T(1:3,1:3,k), 'fro') % this is zero (machine precision)
%     
%     % The translational part differs
%     norm(transls(:,k) - T(1:3,4,k)) %  starts and edns at 0 but may be non-zero in between
% end

%% Plot the interpolated axes in SO(3) x R^3 on the same plot as other axes
figure(1),
for k = 1:num_poses
    T_so3r3 = [Rots(:,:,k), transls(:,k); 0 0 0 1];
    p = T_so3r3 * xyz_homog;
    hold on,
    plot3(p(1,1:2),p(2,1:2),p(3,1:2),'m')  % plot x axis
    plot3(p(1,[1,3]),p(2,[1,3]),p(3,[1,3]),'y') % plot y axis
    plot3(p(1,[1,4]),p(2,[1,4]),p(3,[1,4]),'c') % plot z axis
    axis vis3d, axis equal, grid on
    drawnow
end
hold off
