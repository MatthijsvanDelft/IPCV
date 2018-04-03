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

% Probability gaussian.
sigma = 30; %Trial and error.
h = fspecial('gaussian', [widthSearchArea heightSearchArea], sigma);
normH = h - min(h(:));
h = normH ./ max(normH(:));

% Blob analyser.
minBlobArea = 6;
maxBlobArea = 20;
blobInfo = vision.BlobAnalysis('LabelMatrixOutputPort', true, 'EccentricityOutputPort', true, 'MinimumBlobArea', minBlobArea, 'MaximumBlobArea', maxBlobArea);

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
        
        % Blob analysis.
        [area, centroid, bbox, eccentricity, labeled] = blobInfo.step(dila);
        imshow(label2rgb(labeled));
        hold on;
        plot(widthSearchArea/2, heightSearchArea/2, 'r+');

        temp = eccentricity < minEccentricity;
        circularEnough = eccentricity(eccentricity < minEccentricity);
        circularIndices = find(eccentricity < minEccentricity);
        
        % Calculate probability
        numberBlobs = size(eccentricity, 1);
        blobProb = zeros(numberBlobs,1);
        
        for b = 1:numberBlobs
            blobProb(b,:) = (1-eccentricity(b,:))*h(round(centroid(b,1)), round(centroid(b,2)));                        
        end
        
        thresProb = 0.005;
        [M,I] = max(blobProb(:));
        if (M >= thresProb)
            xBuoy = round(centroid(I,1)) + frameRef(1);
            yBuoy = round(centroid(I,2)) + frameRef(2);                
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