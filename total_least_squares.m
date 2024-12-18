source "./geometry_helpers_2d.m"
source "./total_least_squares_indices.m"
source "./total_least_squares_poses.m"
source "./total_least_squares_projections.m"

# implementation of the boxplus
# applies a perturbation to a set of landmarks and robot poses
# input:
#   XR: the robot poses (3x3xnum_poses: array of homogeneous matrices)
#   XL: the landmark pose (3xnum_landmarks matrix of landmarks)
#   num_poses: number of poses in XR (added for consistency)
#   num_landmarks: number of landmarks in XL (added for consistency)
#   dx: the perturbation vector of appropriate dimensions
#       the poses come first, then the landmarks
# output:
#   XR: the robot poses obtained by applying the perturbation
#   XL: the landmarks obtained by applying the perturbation

function [XR, XL]=boxPlus(XR, XL, num_poses, num_landmarks, dx, pose_dim, landmark_dim)
  for(pose_index=1:num_poses)
    pose_matrix_index=poseMatrixIndex(pose_index, num_poses, num_landmarks);
    dxr=dx(pose_matrix_index:pose_matrix_index+pose_dim-1);
    XR(:,:,pose_index)=v2t(dxr)*XR(:,:,pose_index);
  endfor;
  for(landmark_index=1:num_landmarks)
    landmark_matrix_index=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks);
    dxl=dx(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,:);
    XL(:,landmark_index)+=dxl;
  endfor;
endfunction;


# implementation of the optimization loop with robust kernel
# applies a perturbation to a set of landmarks and robot poses
# input:
#   XR: the initial robot poses (3x3xnum_poses: array of homogeneous matrices)
#   XL: the initial landmark estimates (3xnum_landmarks matrix of landmarks)
#   Z:  the measurements (3xnum_measurements)
#   associations: 2xnum_measurements. 
#                 associations(:,k)=transpose([p_idx,l_idx]) means the kth measurement
#                 refers to an observation made from pose p_idx, that
#                 observed landmark l_idx
#   num_poses: number of poses in XR (added for consistency)
#   num_landmarks: number of landmarks in XL (added for consistency)
#   num_iterations: the number of iterations of least squares
#   damping:      damping factor (in case system not spd)
#   kernel_threshod: robust kernel threshold

# output:
#   XR: the robot poses after optimization
#   XL: the landmarks after optimization
#   chi_stats_{p,r}: array 1:num_iterations, containing evolution of chi2 for projections and poses
#   num_inliers{p,r}: array 1:num_iterations, containing evolution of inliers projections and poses

function [XR, XL, chi_stats_p, num_inliers_p,chi_stats_r, num_inliers_r, H, b] = doTotalLS(XR, XL,
	     Zp, projection_associations,
	     Zr, pose_associations,
	     num_poses,
	     num_landmarks,
	     num_iterations,
	     damping,
	     kernel_threshold,
       kernel_threshold_p,
       pose_dim,
       landmark_dim, K,image_rows, image_cols, rTc)

  chi_stats_p=zeros(1,num_iterations);
  num_inliers_p=zeros(1,num_iterations);
  chi_stats_r=zeros(1,num_iterations);
  num_inliers_r=zeros(1,num_iterations);

  # size of the linear system
  system_size=pose_dim*num_poses+landmark_dim*num_landmarks; 
  for (iteration=1:num_iterations)
    H=zeros(system_size, system_size);
    b=zeros(system_size,1);
   
    if (num_landmarks) 
      [H_projections, b_projections, chi_, num_inliers_] = linearizeProjections(XR, XL, Zp, projection_associations,num_poses, num_landmarks, kernel_threshold_p, pose_dim, landmark_dim,K,image_rows, image_cols, rTc);
      chi_stats_p(iteration)+=chi_;
      num_inliers_p(iteration)=num_inliers_;
    endif;

    [H_poses, b_poses, chi_, num_inliers_] = linearizePoses(XR, XL, Zr, pose_associations,num_poses, num_landmarks, kernel_threshold, pose_dim, landmark_dim);
    chi_stats_r(iteration)+=chi_;
    num_inliers_r(iteration)=num_inliers_;
    
    
    H=H_poses;
    b=b_poses;
    if (num_landmarks) 
       H+=H_projections;
       b+=b_projections;
    endif;
    
    H+=eye(system_size)*damping;
    dx=zeros(system_size,1);

    % we solve the linear system, blocking the first pose
    % this corresponds to "remove" from H and b the locks
    % of the 1st pose, while solving the system
    XL_coords_only = XL(1:3, :); 
    dx(pose_dim+1:end)=-(H(pose_dim+1:end,pose_dim+1:end)\b(pose_dim+1:end,1));
    [XR, XL_coords_only]=boxPlus(XR,XL_coords_only,num_poses, num_landmarks, dx, pose_dim, landmark_dim);
    XL(1:3, :) = XL_coords_only;

  endfor
endfunction
