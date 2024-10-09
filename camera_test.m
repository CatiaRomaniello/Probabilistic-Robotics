source "/home/catia/Probabilistic_Robotics/Probabilistic-Robotics/geometry_helpers_3d.m"
%source "/home/catia/Probabilistic_Robotics/Probabilistic-Robotics/geometry_helpers_2d.m"

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
file_path = '/home/catia/Probabilistic_Robotics/Probabilistic-Robotics/03-PlanarMonocularSLAM/data/world.dat';  
data = dlmread(file_path);

landmark_ids = data(:, 1);       
world_points = data(:, 2:4); 


trajectory_path = '/home/catia/Probabilistic_Robotics/Probabilistic-Robotics/03-PlanarMonocularSLAM/data/trajectoy.dat';
trajectory_data = dlmread(trajectory_path);
pose_ids = trajectory_data(:,1);
odometry_poses = trajectory_data(:, 2:4);  
groundtruth_poses = trajectory_data(:, 5:7);  
%disp(groundtruth_poses);

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
%disp(measurements{2}.gt_pose)
%disp(measurements{51}.image_points)


function point_3D = triangulate(first_obs, second_obs, P1, P2)

    u1 = first_obs(1); 
    v1 = first_obs(2);  
    
    u2 = second_obs(1);  
    v2 = second_obs(2);  

    A = [
        (u1 * P1(3,:) - P1(1,:));  
        (v1 * P1(3,:) - P1(2,:));  
        (u2 * P2(3,:) - P2(1,:));  
        (v2 * P2(3,:) - P2(2,:))  
    ];

    [~, ~, V] = svd(A);
    
    X_homogeneous = V(:,end);
    
    
    point_3D = X_homogeneous(1:3) ./ X_homogeneous(4);
end



unresolved_landmarks = {};
triangulated_points  = {};
for i = 1:length(measurements)
    meas_curr = measurements{i};
    
    odom_curr = meas_curr.odom_pose;
    wTr_curr = v2t([odom_curr(1), odom_curr(2), 0, 0, 0, odom_curr(3)]);
    wTc_curr = wTr_curr * rTc;  
    
    for j = 1:length(meas_curr.actual_ids)
        actual_id = meas_curr.actual_ids(j);
        
        
        is_solved = any(cellfun(@(x) x.actual_id == actual_id, triangulated_points ));
        
        if is_solved
            continue;
        end
        
        idx_unresolved = find(cellfun(@(x) x.actual_id == actual_id, unresolved_landmarks));
        
        if ~isempty(idx_unresolved)
            
            first_obs = unresolved_landmarks{idx_unresolved}.first_observation;
            first_pose = unresolved_landmarks{idx_unresolved}.first_pose;
            
            img_point = meas_curr.image_points(j:j+1, :);
            second_obs = transpose(img_point);
   
            P1 = K * [first_pose(1:3, 1:3), first_pose(1:3, 4)];
            P2 = K * [wTc_curr(1:3, 1:3), wTc_curr(1:3, 4)];
            
            landmark_3D = triangulate(first_obs, second_obs, P1, P2);

            triangulated_points {end + 1} = struct('actual_id', actual_id, 'position', landmark_3D);

            unresolved_landmarks(idx_unresolved) = [];
        else
            
            img_point = meas_curr.image_points(j:j+1, :);
            img_point = transpose(img_point);
            

            unresolved_landmarks{end + 1} = struct( ...
                'actual_id', actual_id, ...
                'first_observation', img_point, ...
                'first_pose', wTc_curr, ...
                'num_frames_unresolved', 0);
        end
    end
    
end


landmark_positions = [];

for i = 1:length(triangulated_points)
    pos_3D = triangulated_points{i}.position;
    
    
    landmark_positions = [landmark_positions; pos_3D];
end

num_landmarks = floor(length(landmark_positions) / 3); 

X = zeros(num_landmarks, 1);
Y = zeros(num_landmarks, 1);
Z = zeros(num_landmarks, 1);

for i = 1:num_landmarks
    idx = (i-1) * 3;  
    X(i) = landmark_positions(idx + 1);  
    Y(i) = landmark_positions(idx + 2);  
    Z(i) = landmark_positions(idx + 3);  
end

figure;
hold on;

plot3(world_points(:,1), world_points(:,2), world_points(:,3), 'bo', 'MarkerSize', 5);  


plot3(X, Y, Z, 'ro', 'MarkerSize', 5);  

xlabel('X');
ylabel('Y');
zlabel('Z');
legend('Ground truth', 'Triangulated landmark');  
grid on;
axis equal;
view(3);  
hold off;
waitfor(gcf); 

