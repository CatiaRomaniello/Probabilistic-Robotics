file_path = '/home/catia/Probabilistic_Robotics/Probabilistic-Robotics/03-PlanarMonocularSLAM/data/world.dat';  
data = dlmread(file_path);


landmark_ids = data(:, 1);
world_points = data(:, 2:4);


figure;
scatter3(world_points(:,1), world_points(:,2), world_points(:,3), 'filled');
hold on;

xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
hold off;

waitfor(gcf); 



