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
widthSearchArea = 100; % In pixels.
heightSearchArea = 100; % In pixels.
FPS = (1/5); % 5 Frames per Second

%% Show the video
% Read in the video and get its width and height.
video = VideoReader('MAH01462.wmv');
videoWidth = video.Width;
videoHeight = video.Height;
videoFPS = video.FrameRate;

% Plot the lens distortion
ut_plot_lens_distortion(cameraParams, [videoHeight videoWidth]);

% Initialise a frame counter.
currentFrame = 0;
% Initialise initial buoy coordinates.
xBuoy = [];
yBuoy = [];

ptThresh = 0.1;

% Create optical flow object using Lucas-Kanade
flowObj = opticalFlowLKDoG( 'NoiseThreshold', 0.002, 'NumFrames', 3,...
                            'ImageFilterSigma', 1.5, ...
                            'GradientFilterSigma', 2);

%Loop through the video.
flowPhaseMag = figure;
videoFigure = figure;
while hasFrame(video)
    tic;
    figure(videoFigure)
    % Video has a new frame, thus increment currentFrame.
    currentFrame = currentFrame + 1;
    frame = readFrame(video, 'native');
    
    % Fix lens distortion.
    [frameUndistorted,~] = undistortImage(frame,cameraParams);
    
    if currentFrame == 1
        % Only extract the coordinates of the buoy in the first frame.
        % GCF is the MATLAB key to the current figure.
        imshow(frame)
        while (size(xBuoy, 1) > 1 || size(yBuoy, 1) > 1) || (isempty(xBuoy) || isempty(yBuoy))
            uiwait(msgbox({'Please pick one point.';...
                           'A shift-, right-, or double-click adds a final point and ends the selection.';...
                           'Pressing Return or Enter ends the selection without adding a final point. ';...
                           'Pressing Backspace or Delete removes the previously selected point.'}, 'Attention!', 'help'));
            [xBuoy, yBuoy] = getpts(gcf);
        end
    end
    
    if currentFrame >=2
        pointsCur = detectFASTFeatures(rgb2gray(frameUndistorted), 'MinContrast', ptThresh);
        pointsPrev = detectFASTFeatures(rgb2gray(framePrev), 'MinContrast', ptThresh);
        
        [featuresCur, pointsCur] = extractFeatures(rgb2gray(frameUndistorted), pointsCur);
        [featuresPrev, pointsPrev] = extractFeatures(rgb2gray(framePrev), pointsPrev);
        
        indexPairs = matchFeatures(featuresPrev, featuresCur);
        pointsPrev = pointsPrev(indexPairs(:,1), :);
        pointsCur = pointsCur(indexPairs(:,2), :);
        
%         showMatchedFeatures(framePrev, frameUndistorted, pointsPrev, pointsCur);
        
        [tform, pointsCurm, pointsPrevm] = estimateGeometricTransform(...
                                                pointsCur, pointsPrev, 'similarity');
        frameUndistortedWarped = imwarp(frameUndistorted, tform, 'OutputView', imref2d(size(frameUndistorted)));
        pointsCurWarped = transformPointsForward(tform, pointsCurm.Location);

        frameCutout = frameUndistortedWarped(yBuoy - 0.5*heightSearchArea : yBuoy + 0.5*heightSearchArea,...
                                             xBuoy - 0.5*widthSearchArea : xBuoy + 0.5*widthSearchArea,...
                                             :);
        
        flow = flowObj.estimateFlow(rgb2gray(frameCutout));
        imshow(rgb2gray(frameCutout))
        hold on
        plot(flow);
        hold off
        figure(flowPhaseMag);
        polarscatter(-flow.Orientation(:), flow.Magnitude(:), '.');
        
%         % Draw the search grid in the image
%         rectangle( 'Position',[xBuoy-0.5*widthSearchArea,...
%                    yBuoy-0.5*heightSearchArea, widthSearchArea,...
%                    heightSearchArea], 'EdgeColor', 'r', 'LineWidth', 3,...
%                    'LineStyle','-', 'Curvature', 0.2)
%         hold off
    end
    framePrev = frameUndistorted;
    T = toc
    
%   Limit to 10 FPS
%     if T < FPS
%         FPS-T
%         pause(0.04)
%     end
    drawnow limitrate
end