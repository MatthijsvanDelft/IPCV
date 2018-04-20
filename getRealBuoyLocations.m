%% Final assignment
% Clear the workspace and close all figures.
clear variables
close all

writeOutputVideo = true;

%% Retrieve calibration information.
% Load the camera parameters.
load('cameraParams.mat');
% Load the estimation errors during the camera calibration.
load('estimationErrors.mat');

% Parameters
widthSearchArea = 100; % In pixels.
heightSearchArea = 100; % In pixels.
keyPointThreshold = 0.1;
distProb_Sigma = 30; % Distance probability (gaussian), Trial and error.
minBlobArea = 6;
maxBlobArea = 20;
adaptThreshSensitivity = 0.3; % Sensitivity of the adapt threshold.
minProb = 0.5; % Minimum probability that detected blob is buoy.
minFlow = 0;
maxFlow = 0.1;
radiusEarth = 6371000; %m
cameraHeight = 2.5; %m
realDistanceHorizon = sqrt(2*radiusEarth*cameraHeight+cameraHeight^2);

%% Read the video and initialise the output video.
% Read in the video and get its width and height.
video = VideoReader('MAH01462.wmv');
videoWidth = video.Width;
videoHeight = video.Height;
videoFPS = video.FrameRate;

if writeOutputVideo
    outVideo = VideoWriter('OutputVideo.avi');
    open(outVideo);
end

%% Initialisation
% Plot the lens distortion for informative purposes.
ut_plot_lens_distortion(cameraParams, [videoHeight videoWidth]);

% Initialise a frame counter.
currentFrame = 227;

% Initialise initial buoy coordinates.
xBuoy = 860.25;
yBuoy = 541.75;
roi(1)=xBuoy; roi(2)=yBuoy;

% Create optical flow object using Lucas-Kanade Derivative of Gaussian
flowObj = opticalFlowLKDoG( 'NoiseThreshold', 0.0039,...
                            'NumFrames', 3,...
                            'ImageFilterSigma', 1.5, ...
                            'GradientFilterSigma', 1);

%Loop through the video.
flowPhaseMag = figure;
videoFigure = figure;

% Probability gaussian.
h = fspecial('gaussian', [widthSearchArea heightSearchArea], distProb_Sigma);
normH = h - min(h(:));
h = normH ./ max(normH(:));

% Blob analyser.
blobInfo = vision.BlobAnalysis('LabelMatrixOutputPort', true,...
                               'EccentricityOutputPort', true,...
                               'MinimumBlobArea', minBlobArea,...
                               'MaximumBlobArea', maxBlobArea,...
                               'ExcludeBorderBlobs', true);
meanDistance = 0;
video.CurrentTime = 227/25;
while hasFrame(video)
    tic;
    %figure(videoFigure)
    % Video has a new frame, thus increment currentFrame.
    currentFrame = currentFrame + 1;
    frame = readFrame(video, 'native');
    prevVect = widthSearchArea;
    closestToCenter = [widthSearchArea/2 heightSearchArea/2];
    
    % Fix lens distortion.
    [frameUndistorted,~] = undistortImage(frame,cameraParams);
    worldMapping = imref2d(round(size(frameUndistorted)*1.4));

    if currentFrame == 228
        % Only extract the coordinates of the buoy in the first frame.
        % GCF is the MATLAB key to the current figure.
        imshow(frame)
%         while (size(xBuoy, 1) > 1 || size(yBuoy, 1) > 1) || (isempty(xBuoy) || isempty(yBuoy))
            uiwait(msgbox({'Please pick one point.';...
                           'A shift-, right-, or double-click adds a final point and ends the selection.';...
                           'Pressing Return or Enter ends the selection without adding a final point. ';...
                           'Pressing Backspace or Delete removes the previously selected point.'}, 'Attention!', 'help'));
            [xBuoy, yBuoy] = getpts(gcf);
            roi(1)=xBuoy; roi(2)=yBuoy;
%         end
        framePrev = frameUndistorted;
    end
    
    if currentFrame >=2
        % Keypoint detection.
        pointsCur = detectFASTFeatures(rgb2gray(frameUndistorted), 'MinContrast', keyPointThreshold);
        pointsPrev = detectFASTFeatures(rgb2gray(framePrev), 'MinContrast', keyPointThreshold);
        
        % Feature extraction.
        [featuresCur, pointsCur] = extractFeatures(rgb2gray(frameUndistorted), pointsCur);
        [featuresPrev, pointsPrev] = extractFeatures(rgb2gray(framePrev), pointsPrev);
        
        % Keypoint matching.
        indexPairs = matchFeatures(featuresPrev, featuresCur);
        pointsPrev = pointsPrev(indexPairs(:,1), :);
        pointsCur = pointsCur(indexPairs(:,2), :);

%         showMatchedFeatures(framePrev, frameUndistorted, pointsPrev, pointsCur);
        
        % Stabilization transformation.
        [tform, pointsCurm, pointsPrevm] = estimateGeometricTransform(pointsCur, pointsPrev, 'similarity');
        [frameUndistortedWarped,~] = imwarp(frameUndistorted, tform, 'OutputView', worldMapping);

        % Memory for camera stabilisation    
        framePrev = frameUndistortedWarped;   
        
        % ROI.
        frameRef = [roi(1) - 0.5*widthSearchArea roi(2) - 0.5*heightSearchArea];
        searchArea = frameUndistortedWarped(roi(2) - 0.5*heightSearchArea : roi(2) + 0.5*heightSearchArea,...
                                             roi(1) - 0.5*widthSearchArea : roi(1) + 0.5*widthSearchArea,...
                                             :);
        imshow(searchArea);
        videoFigure.WindowState = 'maximized';
        temp1 = []; temp2 = [];
       while (size(temp1, 1) > 1 || size(temp2, 1) > 1) || (isempty(temp1) || isempty(temp2))
            [temp1, temp2] = getpts(gcf);

            if isempty(temp1) || isempty(temp2)
                xBuoy(currentFrame) = -1; yBuoy(currentFrame) = -1;
                temp1 = -1; temp2 = -1;
            else (size(temp1, 1) == 1 || size(temp2, 1) == 1)
                xBuoy(currentFrame) = round(temp1) + frameRef(1);
                yBuoy(currentFrame) = round(temp2) + frameRef(2); 
                roi(1)=xBuoy(currentFrame); roi(2)=yBuoy(currentFrame);
            end
       end
    end
    
    proc_time = toc
    drawnow limitrate
end

% Close the output video to actually make the data available outside
% matlab.
if writeOutputVideo
    close(outVideo);
end