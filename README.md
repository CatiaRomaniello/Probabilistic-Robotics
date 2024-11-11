# Probabilistic-Robotics
Planar monocular SLAM project.

The following project aims to build a map where the robot moves while localizing itself (SLAM).
The differential drive robot is equipped with a monocular camera.
Given the wheeled odometry, the stream of point projections, and the camera parameters as input, the objective is to estimate the robot's trajectory and the 3D landmarks' pose. These results are evaluated by computing the Root Mean Square Error (RMSE) between the estimates and the corresponding ground truth values.
## How to compile and run
Clone the repository
```
git clone https://github.com/CatiaRomaniello/Probabilistic-Robotics.git
```
and select the folder 

```
cd Probabilistic-Robotics/
```

finally, to run it 

```
cd octave TotalLeastSquarenew.m
```
## Results

By running the code, the terminal will show the RMSE before and after the optimization for both the robot pose(rotational and translational error) and the positions of the landmarks.
After a few minutes, it will show some plots comparing the results before (the initial guess and the ground truth) and after (the optimized results and the ground truth) the optimization.

Below is a plot of the robot's trajectory considering the initial guess, the optimized result, and the ground truth.
For more details and results, please refer to the **PR_report.pdf** file.

![Screenshot from 2024-11-06 20-35-24](https://github.com/user-attachments/assets/b0f50e24-01c8-4dda-9185-6723bd02ac36)


![Screenshot from 2024-11-11 20-42-29](https://github.com/user-attachments/assets/cff9d582-1fb7-4b59-8804-eec8af439cf7)

Once the code is started, as soon as the first plot is displayed, press enter to view all others. 
If you press enter after Figure 5 is generated, it will close all the outputs.