function posevecs = PoseMats2PoseVecs(PoseMats)

num_poses = size(PoseMats,3);
posevecs = nan*ones(6,num_poses);
for k = 1:num_poses
    posevecs(:,k) = vee_se3( logm( PoseMats(:,:,k) ));
end
