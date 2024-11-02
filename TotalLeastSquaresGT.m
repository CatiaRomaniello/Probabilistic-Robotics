source "./total_least_squares.m"

pose_dim = 3;  
landmark_dim = 3; 
file_path = '/home/catia/Probabilistic_Robotics/Probabilistic-Robotics/03-PlanarMonocularSLAM/data/world.dat';  
data = dlmread(file_path);

landmark_ids = data(:, 1);  
XL_true = transpose(data(:, 2:4));  
num_landmarks = size(XL_true, 2);  
file_path = '/home/catia/Probabilistic_Robotics/Probabilistic-Robotics/03-PlanarMonocularSLAM/data/camera.dat';  
data = dlmread(file_path);
K = data(2:4,1:3);
rTc = data(6:9,1:4);
z_n= data(10,2);
z_f = data(11,2);
image_cols = data(12,2);
image_rows = data(13,2);
folder_path = '/home/catia/Probabilistic_Robotics/Probabilistic-Robotics/03-PlanarMonocularSLAM/data/';
file_base_name = 'meas-';  
file_extension = '.dat';   

measurements = cell(1, 200);  

for i = 0:199
    file_number = sprintf('%05d', i);  
    meas_path = fullfile(folder_path, [file_base_name, file_number, file_extension]);
    
    fid = fopen(meas_path, 'r');  
    
    if fid == -1
        warning(['Impossibile aprire il file: ', meas_path]);
        continue;  
    end
    
    line = fgetl(fid);
    seq_numb = sscanf(line, 'seq: %d');
    
    line = fgetl(fid);
    gt_pose = sscanf(line, 'gt_pose: %f %f %f');
    
    line = fgetl(fid);
    odom_pose = sscanf(line, 'odom_pose: %f %f %f');
    
    point_ids = [];
    actual_ids = [];
    image_points = [];
    
    while ~feof(fid)
        line = fgetl(fid);
        if strncmp(line, 'point', 5)
            data = sscanf(line, 'point %d %d %f %f');
            point_ids = [point_ids; data(1)];
            actual_ids = [actual_ids; data(2)];
            image_points = [image_points; data(3:4)];
        end
    end
    
    fclose(fid);
    
    s = struct('seq_numb', seq_numb, 'gt_pose', gt_pose, 'odom_pose', odom_pose, ...
               'point_ids', point_ids, 'actual_ids', actual_ids, 'image_points', image_points);
    
    measurements{i + 1} = s;  
    
end

trajectory_path = '/home/catia/Probabilistic_Robotics/Probabilistic-Robotics/03-PlanarMonocularSLAM/data/trajectoy.dat';
trajectory_data = dlmread(trajectory_path);

pose_ids = trajectory_data(:, 1);
odometry_poses = trajectory_data(:, 2:4);  
groundtruth_poses = trajectory_data(:, 5:7); 

num_poses = size(pose_ids, 1);
XR_true = zeros(3, 3, num_poses);  

for i = 1:num_poses
    XR_true(:,:,i) = v2t([groundtruth_poses(i,1), groundtruth_poses(i,2), groundtruth_poses(i,3)]); 
end
for i = 1:num_poses
    XR_guess(:,:,i) = v2t([odometry_poses(i,1), odometry_poses(i,2), odometry_poses(i,3)]);  
end




######################################## PROJECTION MEASUREMENTS ########################################


measurement_num = 1;  
num_projection_measurements = 1;  

for i = 1:length(measurements)
    meas = measurements{i}; 
    
    if isempty(meas)  
        continue;
    end
    
    pose_num = meas.seq_numb + 1; 
    point_ids = meas.point_ids; 
    actual_ids = meas.actual_ids; 
    image_points = meas.image_points;  
    
    k = 1;
    for j = 1:length(point_ids)
        landmark_num = actual_ids(j)+1;  
        z_img = image_points(k:k+1, :); 
        projection_associations(:, measurement_num) = [pose_num; landmark_num];
        Zp(:, measurement_num) = transpose(z_img);
        measurement_num++;
        k+=2;
    end
end




######################################## POSE MEASUREMENTS ########################################



num_pose_measurements = num_poses - 1;
Zr = zeros(3, 3, num_pose_measurements);
pose_associations = zeros(2, num_pose_measurements);

for pose_num = 1:num_pose_measurements
    Xi = v2t([groundtruth_poses(pose_num, 1), groundtruth_poses(pose_num, 2), groundtruth_poses(pose_num, 3)]);
    Xj = v2t([groundtruth_poses(pose_num + 1, 1), groundtruth_poses(pose_num + 1, 2), groundtruth_poses(pose_num + 1, 3)]);
    pose_associations(:, pose_num) = transpose([pose_num, pose_num + 1]);
    Zr(:, :, pose_num) = inv(Xi)*Xj;
end
#for pose_num = 1:num_pose_measurements
#    Xi = v2t([odometry_poses(pose_num, 1), odometry_poses(pose_num, 2), odometry_poses(pose_num, 3)]);
#    Xj = v2t([odometry_poses(pose_num + 1, 1), odometry_poses(pose_num + 1, 2), odometry_poses(pose_num + 1, 3)]);
#    pose_associations(:, pose_num) = transpose([pose_num, pose_num + 1]);
#    Zr(:, :, pose_num) = inv(Xi) * Xj;
#end

############################## GROUNDTRUTH INITIAL GUESS ################################## 

XR_guess=XR_true;
XL_guess=XL_true;



############################## CALL SOLVER  ################################## 

# uncomment the following to suppress pose-landmark-projection measurements
#num_landmarks=0;
#Zp=zeros(3,0);

# uncomment the following to suppress pose-pose measurements
%Zr=zeros(3,3,0);

damping=1;
kernel_threshold=1;
kernel_threshold_p=1000;
num_iterations=30;
[XR, XL, chi_stats_p, num_inliers_p, chi_stats_r, num_inliers_r, H, b]=doTotalLS(XR_guess, XL_guess, 
												      Zp, projection_associations, 
												      Zr, pose_associations, 
												      num_poses, 
												      num_landmarks, 
												      num_iterations, 
												      damping, 
												      kernel_threshold,
                                                      kernel_threshold_p,
                                                      pose_dim, 
                                                      landmark_dim,K,image_rows, image_cols, rTc);
												      





rotation_errors = [];
translation_errors = [];

for i = 2:num_poses

    rel_T_optimized = inv(XR(:, :, i-1)) * XR(:, :, i);
    rel_T_gt = inv(XR_guess(:, :, i-1)) * XR_guess(:, :, i);
    error_T = inv(rel_T_optimized) * rel_T_gt;
    rot_error = atan2(error_T(2, 1), error_T(1, 1));
    rotation_errors = [rotation_errors; rot_error];
    trans_error = norm(error_T(1:2, 3));  % Errore di traslazione (solo X e Y)
    translation_errors = [translation_errors; trans_error];
end

rmse_rotation = sqrt(sum(rotation_errors .^ 2) / length(rotation_errors));
rmse_translation = sqrt(sum(translation_errors .^ 2) / length(translation_errors));

disp(rmse_rotation);
disp(rmse_translation);

errors_guess = sqrt(sum((XL_guess - XL_true).^2, 1)); 
rmse_landmarks_guess = sqrt(mean(errors_guess.^2));    
disp('RMSE of landmarks before optimization:');
disp(rmse_landmarks_guess);


errors_optimized = sqrt(sum((XL - XL_true).^2, 1));    
rmse_landmarks_optimized = sqrt(mean(errors_optimized.^2)); 
disp('RMSE of landmarks after optimization:');
disp(rmse_landmarks_optimized);
     
												      
# Plot State
XR_vec = zeros(3, num_poses);
XR_true_vec = zeros(3, num_poses);
XR_guess_vec = zeros(3, num_poses);

for i = 1:num_poses
    XR_vec(:, i) = t2v(XR(:, :, i));
    XR_true_vec(:, i) = t2v(XR_true(:, :, i));
    XR_guess_vec(:, i) = t2v(XR_guess(:, :, i));
end
figure(1);
hold on;
grid;

subplot(2,2,1);
title("Landmark Initial Guess");
plot(XL_true(1,:),XL_true(2,:),'b*',"linewidth",2);
hold on;
plot(XL_guess(1,:),XL_guess(2,:),'ro',"linewidth",2);
legend("Landmark True", "Guess");grid;


subplot(2,2,2);
title("Landmark After Optimization");
plot(XL_true(1,:),XL_true(2,:),'b*',"linewidth",2);
hold on;
plot(XL(1,:),XL(2,:),'ro',"linewidth",2);
legend("Landmark True", "Guess");grid;


subplot(2,2,3);
title("Poses Initial Guess");
plot(XR_true_vec(1,:),XR_true_vec(2,:),'b*',"linewidth",2);
hold on;
plot(XR_guess_vec(1,:),XR_guess_vec(2,:),'ro',"linewidth",2);
legend("Poses True", "Guess");grid;


subplot(2,2,4);
title("Poses After Optimization");
plot(XR_true_vec(1,:),XR_true_vec(2,:),'b*',"linewidth",2);
hold on;
plot(XR_vec(1,:),XR_vec(2,:),'ro',"linewidth",2);
legend("Poses True", "Guess"); grid;
pause(); 

figure(2);
hold on;
grid;
title("chi evolution");

subplot(3,2,1);
plot(chi_stats_r, 'r-', "linewidth", 2);
legend("Chi Poses"); grid; xlabel("iterations");
subplot(3,2,2);
plot(num_inliers_r, 'b-', "linewidth", 2);
legend("#inliers"); grid; xlabel("iterations");

subplot(3,2,3);
plot(chi_stats_p, 'r-', "linewidth", 2);
legend("Chi Proj"); grid; xlabel("iterations");
subplot(3,2,4);
plot(num_inliers_p, 'b-', "linewidth", 2);
legend("#inliers");grid; xlabel("iterations");
pause(); 
figure(3);
title("H matrix");
H_ =  H./H;                      # NaN and 1 element
H_(isnan(H_))=0;                 # Nan to Zero
H_ = abs(ones(size(H_)) - H_);   # switch zero and one
H_ = flipud(H_);                 # switch rows
colormap(gray(64));
hold on;
image([0.5, size(H_,2)-0.5], [0.5, size(H_,1)-0.5], H_*64);
hold off;
pause(); 