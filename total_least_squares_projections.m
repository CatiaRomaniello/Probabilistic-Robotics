# Assembly of the projection problem
source "./geometry_helpers_3d.m"
source "./geometry_helpers_2d.m"
source "./total_least_squares_indices.m"
# error and jacobian of a measured landmark
# input:
#   Xr: the robot pose in world frame (3x3 homogeneous matrix) ->then transformed in 4x4 homogenepus matrix since rTc is 4x4
#   Xl: the landmark pose (3x1 vector, 3d pose in world frame)
#   z:  projection of the landmark on the image plane
# output:
#   e: 2x1 is the difference between prediction and measurement
#   Jr: 2x3 derivative w.r.t a the error and a perturbation on the
#       pose
#   Jl: 2x3 derivative w.r.t a the error and a perturbation on the
#       landmark
#   is_valid: true if projection ok
# K = camera matrix
# rTc = camera frame wrt to robot frame

function [is_valid, e,Jr,Jl]=projectionErrorAndJacobian(Xr,Xl,z,K, image_rows, image_cols, rTc)
  is_valid=false;
  e=[0;0];
  Jr=zeros(2,3);
  Jl=zeros(2,3);
  Xr_3d = eye(4);
  Xr_3d(1:2,1:2) = Xr(1:2,1:2);
  Xr_3d(1:2,4) = Xr(1:2,3);

  wTc = Xr_3d*rTc; #camera frame wrt to world frame 
  iR = transpose(wTc(1:3,1:3));
  it = -iR*wTc(1:3,4);

  pw=iR*Xl+it; 

  if (pw(3)<0)
     return;
  endif

  Jwr_t1=zeros(3,6);
  Jwr_t1(1:3,1:3)=-iR;
  Jwr_t1(:,4:6)=iR*skew(Xl);

  Jwr = zeros(3,3);
  Jwr(1:3,1:2) = Jwr_t1(1:3,1:2);
  Jwr(1:3,3) = Jwr_t1(1:3,6);
  Jwl=iR;
  
  p_cam=K*pw;
  iz=1./p_cam(3);
  z_hat=p_cam(1:2)*iz;
  if (z_hat(1)<0 || 
      z_hat(1)>image_cols ||
      z_hat(2)<0 || 
      z_hat(2)>image_rows)
    return;
  endif;

  iz2=iz*iz;
  Jp=[iz, 0, -p_cam(1)*iz2;
      0, iz, -p_cam(2)*iz2];
  
  e=z_hat-z;
  Jr=Jp*K*Jwr;
  Jl=Jp*K*Jwl;
  is_valid=true;
  
endfunction;


#linearizes the robot-landmark measurements
#   XR: the initial robot poses (3x3xnum_poses: array of homogeneous matrices)
#   XL: the initial landmark estimates (3xnum_landmarks matrix of landmarks)
#   Z:  the measurements (2xnum_measurements)
#   associations: 2xnum_measurements. 
#                 associations(:,k)=transpose([p_idx,l_idx]) means the kth measurement
#                 refers to an observation made from pose p_idx, that
#                 observed landmark l_idx
#   num_poses: number of poses in XR (added for consistency)
#   num_landmarks: number of landmarks in XL (added for consistency)
#   kernel_threshod: robust kernel threshold
# output:
#   XR: the robot poses after optimization
#   XL: the landmarks after optimization
#   chi_stats: array 1:num_iterations, containing evolution of chi2
#   num_inliers: array 1:num_iterations, containing evolution of inliers

function [H,b, chi_tot, num_inliers]=linearizeProjections(XR, XL, Zl, associations,num_poses, num_landmarks, kernel_threshold, pose_dim, landmark_dim,K,image_rows, image_cols, rTc)

  system_size=pose_dim*num_poses+landmark_dim*num_landmarks; 
  H=zeros(system_size, system_size);
  b=zeros(system_size,1);
  chi_tot=0;
  num_inliers=0;
  for (measurement_num=1:size(Zl,2))
    pose_index=associations(1,measurement_num);
    actual_id=associations(2,measurement_num);
    landmark_index = find(XL(4, :) == actual_id, 1);
    if isempty(landmark_index) 
    continue; % Salta questo landmark se non è valido
    end

    z=Zl(:,measurement_num);
    Xr=XR(:,:,pose_index);
    Xl = XL(1:3, landmark_index); 
    [is_valid, e,Jr,Jl] = projectionErrorAndJacobian(Xr, Xl, z, K,image_rows, image_cols, rTc);
    if (! is_valid)
       continue;
    endif;
    chi=transpose(e)*e;
    w = 1;
    if (chi>kernel_threshold)
      e*=sqrt(kernel_threshold/chi);
      w = 1/kernel_threshold;
      chi=kernel_threshold;
    else
      num_inliers++;
    endif;
    chi_tot+=chi;

    omega = w*eye(2);
    Hrr = transpose(Jr) * omega * Jr;
    Hrl = transpose(Jr) * omega * Jl;
    Hlr = transpose(Jl) * omega * Jr;
    Hll = transpose(Jl) * omega * Jl;
    br = transpose(Jr) * omega * e;
    bl = transpose(Jl) * omega * e;

    pose_matrix_index=poseMatrixIndex(pose_index, num_poses, num_landmarks);
    
    landmark_matrix_index=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks);



    H(pose_matrix_index:pose_matrix_index+pose_dim-1,
      pose_matrix_index:pose_matrix_index+pose_dim-1)+=Hrr;

    H(pose_matrix_index:pose_matrix_index+pose_dim-1,
      landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=Hrl;

    H(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,
      landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=Hll;

    H(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,
      pose_matrix_index:pose_matrix_index+pose_dim-1)+=Hlr;

    b(pose_matrix_index:pose_matrix_index+pose_dim-1)+=br;
    b(landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=bl;
  endfor
endfunction
