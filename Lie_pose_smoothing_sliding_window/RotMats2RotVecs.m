function rotvecs = RotMats2RotVecs(RotMats)

num_poses = size(RotMats,3);
rotvecs = nan*ones(3,num_poses);
for k = 1:num_poses
    rotvecs(:,k) = vee_so3( logm( RotMats(:,:,k) ));
end
