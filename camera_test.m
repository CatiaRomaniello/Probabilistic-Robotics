file_path = '/home/catia/Probabilistic_Robotics/Probabilistic-Robotics/03-PlanarMonocularSLAM/data/camera.dat';  
data = dlmread(file_path);
K = data(2:4,1:3);
rTc = data(6:9,1:4);
z_n= data(10,2);
z_f = data(11,2);
wid = data(12,2);
hg = data(13,2);

%disp(K);
%disp(rTc);
%disp(z_n);
%disp(z_f);
%disp(wid);
%disp(hg);

trajectory_path = '/home/catia/Probabilistic_Robotics/Probabilistic-Robotics/03-PlanarMonocularSLAM/data/trajectoy.dat';
trajectory_data = dlmread(trajectory_path);
pose_ids = trajectory_data(:,1);
odometry_poses = trajectory_data(:, 2:4);  
groundtruth_poses = trajectory_data(:, 5:7);  
disp(groundtruth_poses);

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



%disp(measurements{51}.seq_numb)
disp(measurements{2}.gt_pose)
%disp(measurements{51}.odom_pose)


