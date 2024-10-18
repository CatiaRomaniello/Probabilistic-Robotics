source "/home/catia/Probabilistic_Robotics/Probabilistic-Robotics/geometry_helpers_3d.m"

%(minimal) size of pose and landmarks
global pose_dim=6;
global landmark_dim=3;


function v_idx=poseMatrixIndex(pose_index, num_poses, num_landmarks)
  global pose_dim;
  global landmark_dim;

  if (pose_index>num_poses)
    v_idx=-1;
    return;
  endif;
  v_idx=1+(pose_index-1)*pose_dim;
endfunction;

function v_idx=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks)
  global pose_dim;
  global landmark_dim;
  if (landmark_index>num_landmarks)
    v_idx=-1;
    return;
  endif;
  v_idx=1 + (num_poses)*pose_dim + (landmark_index-1) * landmark_dim;
endfunction;

function [is_behind, z_hat] = project(Xr, Xl, K)
  p_camera_frame = Xr(1:3,1:3) * Xl + Xr(1:3,4);
  p_projected = K * p_camera_frame;
  is_behind = true;
  z_hat = [-1;-1];
  if(p_camera_frame(3) > 0)
    z_hat = p_projected(1:2)/p_projected(3);
    is_behind = false;
  endif
endfunction


function [associations, Z, seen_landmarks] = computeMeasurementAndAssociations(XR, XL, K)
  num_poses=size(XR,3);
  num_landmarks=size(XL,2);
  associations = zeros(2,0);
  Z = zeros(2,0);
  seen_landmarks = zeros(num_landmarks,1);
  for (pose_num=1:num_poses)
      Xr=XR(:,:,pose_num);
      R = Xr(1:3,1:3);
      t = Xr(1:3,4);
      for (landmark_num=1:num_landmarks)
        Xl=XL(:,landmark_num);
        p = (R*Xl + t);
        if(p(3) < 0)
          continue;
        endif
        z_proj = K * p;
        associations(:,end+1)=transpose([pose_num,landmark_num]);
        # binary array that take track of how many times each landmark have been observed
        seen_landmarks(landmark_num) += 1;
        Z(:,end+1)= z_proj(1:2)/z_proj(3);
      endfor;
  endfor
endfunction


function [is_behind,e,Jr,Jl]=errorAndJacobian(Xr, Xl, z, K)
   R=Xr(1:3,1:3);
   t=Xr(1:3,4);
   [is_behind,z_hat]=project(Xr,Xl,K); #prediction
   e = zeros(2,1);
   Jl = zeros(2,3);
   Jr = zeros(2,6);
   if (is_behind)
    return
   endif
   p_camera_frame = R * Xl + t;
   p_projected = K * p_camera_frame;
   e=z_hat-z;
   Jr=zeros(2,6);
   Jl=zeros(2,3);
   # TODO: fill the jacobians
   inverse_z = 1/p_projected(3);
   inverse_square_z = inverse_z * inverse_z;
   J_proj = [inverse_z, 0, -p_projected(1) * inverse_square_z; 
             0, inverse_z, -p_projected(2) * inverse_square_z];

   %Jicp = J_icp(p_camera_frame);
   Jicp = [eye(3), -skew(p_camera_frame)];
   Jr = J_proj * K * Jicp;
   Jl = J_proj * K * R;
endfunction;

function [XR, XL]=boxPlus(XR, XL, num_poses, num_landmarks, dx)
  global pose_dim;
  global landmark_dim;
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
#   XR: the initial robot poses (4x4xnum_poses: array of homogeneous matrices)
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
#   chi_stats: array 1:num_iterations, containing evolution of chi2
#   num_inliers: array 1:num_iterations, containing evolution of inliers
function [XR, XL, chi_stats, num_inliers]=doBundleAdjustment(XR, XL, Z, K,
							associations, 
							num_poses, 
							num_landmarks, 
							num_iterations, 
							damping, 
							kernel_threshold)
  global pose_dim;
  global landmark_dim;
  

  chi_stats=zeros(1,num_iterations);
  num_inliers=zeros(1,num_iterations);
  # size of the linear system
  system_size=pose_dim*num_poses+landmark_dim*num_landmarks; 
  for (iteration=1:num_iterations)
    H=zeros(system_size, system_size);
    b=zeros(system_size,1);
    chi_stats(iteration)=0;
    disp(['--- Iterazione: ', num2str(iteration), ' ---']);
    for (measurement_num=1:size(Z,2))
      pose_index=associations(1,measurement_num);
      landmark_index=associations(2,measurement_num);
      
      z=Z(:,measurement_num);
      Xr=XR(:,:,pose_index);
      Xl=XL(:,landmark_index);
      if landmark_index <= 0 || landmark_index > num_landmarks
                error('landmark_index non valido: %d', landmark_index);
      endif
            
      [is_behind, e,Jr,Jl] = errorAndJacobian(Xr, Xl, z, K);
      if (is_behind)
        continue;
      endif
      chi=transpose(e)*e;
      # apply the robust kernel as a weight to the linear system
      w = 1;
      if (chi>kernel_threshold)
      	w = 1/kernel_threshold;
      	chi=kernel_threshold;
      else
      	num_inliers(iteration)++;
      endif;
      chi_stats(iteration)+=chi;
    
      omega = w * eye(2);
      Hrr = transpose(Jr) * omega * Jr;
      Hrl = transpose(Jr) * omega * Jl;
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
	pose_matrix_index:pose_matrix_index+pose_dim-1)+=transpose(Hrl);

      b(pose_matrix_index:pose_matrix_index+pose_dim-1)+=br;
      b(landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=bl;

    endfor
    size_H=size(H,1)
    rank_H=rank(H)
    H+=eye(system_size)*damping;
    dx=zeros(system_size,1);
    
    % we solve the linear system, blocking the first pose
    % this corresponds to "remove" from H and b the locks
    % of the 1st pose, while solving the system

    dx(pose_dim+1:end)=-(H(pose_dim+1:end,pose_dim+1:end)\b(pose_dim+1:end,1));
    [XR, XL]=boxPlus(XR,XL,num_poses, num_landmarks, dx);
  endfor

  % visualize the sparsity pattern of the approximate hessian
  visual_H = ones(system_size,system_size);
  for r=1:system_size
    for c=1:system_size
      if(abs(H(r,c)) > 0)
        visual_H(r,c) = 0;
      endif
    endfor
  endfor
  imshow(visual_H)
  waitfor(gcf); 
endfunction

function i = plotState(XL, XL_guess, XL_gt)
#plot landmarks
hold on;
plot3(XL(1,:),XL(2,:),XL(3,:),'b*',"linewidth",2);
hold on;
plot3(XL_guess(1,:),XL_guess(2,:),XL_guess(3,:),'ro',"linewidth",2);
hold on;
plot3(XL_gt(1,:),XL_gt(2,:),XL_gt(3,:),'g*',"linewidth",2);
hold on;
legend("estimate","initial guess","ground truth")
i = 1;
endfunction
