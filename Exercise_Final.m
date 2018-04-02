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
minEccentricity = 0.3; % For perfect circle, eccentricity is 0.

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

blobInfo = vision.BlobAnalysis('LabelMatrixOutputPort', true, 'EccentricityOutputPort', true);
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
        
        frameRef = [xBuoy - 0.5*widthSearchArea yBuoy - 0.5*heightSearchArea];
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
        title('Original image with camera stabilisation');
        
        subplot(2,2,4)
        %imshow(frame)
        %imshow(cutoutFilter);
        Thres = adaptthresh(rgb2gray(cutoutFilter), 0.20);
        bin = imclearborder(imbinarize(rgb2gray(cutoutFilter), Thres));
        erod = imerode(bin, strel('disk', 1));
        dila = imdilate(erod, strel('disk', 1));
        %imshow((dila));
        [area, centroid, bbox, eccentricity, labeled] = blobInfo.step(dila);
        imshow(label2rgb(labeled));
        hold on;
        plot(widthSearchArea/2, heightSearchArea/2, 'r+');
        temp = eccentricity < minEccentricity;
        for i = 1:size(eccentricity, 1)
            if temp(i) == 1
                rectangle('Position', [bbox(i, 1), bbox(i, 2), ...
                          bbox(i, 3), bbox(i, 4)], ...
                          'EdgeColor', 'r',...
                          'LineStyle','-');
                      
                vect = norm( centroid(i,:) - [widthSearchArea/2 heightSearchArea/2]);
                if vect < prevVect
                    prevVect = vect;
                    closestToCenter = centroid(i,:);
                end
            end
        end
        if sum(eccentricity) > 0
            line([widthSearchArea/2 closestToCenter(1)], [heightSearchArea/2 closestToCenter(2)]); % syntax is [x1, x2], [y1, y2]
            xBuoy = frameRef(1) + closestToCenter(1); yBuoy = frameRef(2) + closestToCenter(2);
        end
        hold off;
        title('Thresholded, eroded, dilated and labeled searchgrid');
        
        subplot(2,2,3)
        %polarscatter(-flow.Orientation(:), flow.Magnitude(:), '.');
        imshow(cutoutFilter);
        
        title('Zoomed searchgrid');
        
        subplot(2,2,2)
        imshow((dila));%imclearborder(imbinarize(rgb2gray(cutoutFilter), Thres)));
        %imshow(frameUndistortedWarped)
        %figure(flowPhaseMag);        
        
        title('Thresholded, eroded and dilated searchgrid');
        
    end
    framePrev = frameUndistorted;
    T = toc
    
    drawnow limitrate
end