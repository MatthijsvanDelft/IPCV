%% Final assignment
% Clear the workspace and close all figures.
clear variables
close all

%% Retrieve calibration information.
% Load the camera parameters.
load('cameraParams.mat');
% Load the estimation errors during the camera calibration.
load('estimationErrors.mat');

%% Show the video
video = VideoReader('MAH01462.MP4');
videoWidth = video.Width;
videoHeight = video.Height;

while hasFrame(video)
    frame = readFrame(video, 'native');
    imshow(frame)
end

%getpts();