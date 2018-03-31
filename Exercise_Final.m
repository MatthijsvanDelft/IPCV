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
heightSearchArea = 200; % In pixels.
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
flowObj = opticalFlowLKDoG( 'NoiseThreshold', 0.0012, 'NumFrames', 3,...
                            'ImageFilterSigma', 3.5, ...
                            'GradientFilterSigma', 4.5);
% flowObj = opticalFlowLK( 'NoiseThreshold', 0.0144);
%Loop through the video.
flowPhaseMag = figure;
videoFigure = figure;
while hasFrame(video)
    tic;
    %figure(videoFigure)
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
        %flow = flowObj.estimateFlow(rgb2gray(frameCutout));
        cutoutFilter = frameCutout; %imgaussfilt(frameCutout, 1);
        
        subplot(2,2,1)
        imshow(frameUndistortedWarped);
        hold on
        drawSearchGrid(xBuoy, yBuoy, widthSearchArea, heightSearchArea);
        hold off
        %rotate(mesh(im2double(rgb2gray(cutoutFilter))), [0 0 1], 90)
        %zlim([0 1])
        %plot(flow);
        %showMatchedFeatures(framePrev, frameUndistorted, pointsPrev, pointsCur);
        subplot(2,2,2)
        %imshow(frame)
        %imshow(cutoutFilter);
        Thres = adaptthresh(rgb2gray(cutoutFilter), 0.20);
        bin = imclearborder(imbinarize(rgb2gray(cutoutFilter), Thres));
        erod = imerode(bin, strel('disk', 1));
        dila = imdilate(erod, strel('disk', 1));
        imshow((dila));
        
        subplot(2,2,3)
        %polarscatter(-flow.Orientation(:), flow.Magnitude(:), '.');
        imshow(cutoutFilter);
        
        subplot(2,2,4)
        Thres = adaptthresh(rgb2gray(cutoutFilter), 0.20);
        imshow(imclearborder(bin));%imclearborder(imbinarize(rgb2gray(cutoutFilter), Thres)));
        %imshow(frameUndistortedWarped)
        %figure(flowPhaseMag);        
        
        
    end
    framePrev = frameUndistorted;
    T = toc
    
    drawnow limitrate
end