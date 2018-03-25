%% Final assignment
% Clear the workspace and close all figures.
clear variables
close all

%% Retrieve calibration information.
% Load the camera parameters.
load('cameraParams.mat');
% Load the estimation errors during the camera calibration.
load('estimationErrors.mat');

% Parameters
widthSearchArea = 200; % In pixels.
heightSearchArea = 100; % In pixels.

%% Show the video
% Read in the video and get its width and height.
video = VideoReader('MAH01462.MP4');
videoWidth = video.Width;
videoHeight = video.Height;
videoFPS = video.FrameRate;

% Plot the lens distortion
ut_plot_lens_distortion(cameraParams, [videoHeight videoWidth]);
figure;

% Initialise a frame counter.
currentFrame = 0;
% Initialise initial buoy coordinates.
xBuoyInitial = [];
yBuoyInitial = [];

%Loop through the video.
while hasFrame(video)
    tic;
    
    % Video has a new frame, thus increment currentFrame.
    currentFrame = currentFrame + 1;
    frame = readFrame(video, 'native');
    imshow(frame)
    hold on;
    
    if currentFrame == 1
        % Only extract the coordinates of the buoy in the first frame.
        % GCF is the MATLAB key to the current figure.
        while (size(xBuoyInitial, 1) > 1 || size(yBuoyInitial, 1) > 1) || (isempty(xBuoyInitial) || isempty(yBuoyInitial))
            uiwait(msgbox({'Please pick one point.';...
                           'A shift-, right-, or double-click adds a final point and ends the selection.';...
                           'Pressing Return or Enter ends the selection without adding a final point. ';...
                           'Pressing Backspace or Delete removes the previously selected point.'}, 'Attention!', 'help'));
            [xBuoyInitial, yBuoyInitial] = getpts(gcf);
        end
    end
    % Draw the search grid in the image
    rectangle('Position',[xBuoyInitial-0.5*widthSearchArea, yBuoyInitial-0.5*heightSearchArea, widthSearchArea, heightSearchArea],...
	'EdgeColor', 'r', 'LineWidth', 3, 'LineStyle','-', 'Curvature', 0.2)
    hold off;
    T = toc;
    
    % Limit to video's fps
    if T < 1/videoFPS
        pause((1/videoFPS)-T)
    end
end