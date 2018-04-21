%% Final assignment
% Clear the workspace and close all figures.
clear variables
close all

writeOutputVideo = true;

%% Retrieve calibration information.
% Load the camera parameters.
load('data/cameraParams.mat');
% Load the estimation errors during the camera calibration.
load('data/estimationErrors.mat');

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
    outVideo = VideoWriter('OutputVideo.mp4', 'MPEG-4');
    open(outVideo);
end

%% Initialisation
% Plot the lens distortion for informative purposes.
ut_plot_lens_distortion(cameraParams, [videoHeight videoWidth]);

% Initialise a frame counter.
currentFrame = 0;

% Initialise initial buoy coordinates.
xBuoy = [];
yBuoy = [];

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

while hasFrame(video)
    tic;
    %figure(videoFigure)
    % Video has a new frame, thus increment currentFrame.
    currentFrame = currentFrame + 1;
    frame = readFrame(video, 'native');
    
    % Fix lens distortion.
    [frameUndistorted,~] = undistortImage(frame,cameraParams);
    worldMapping = imref2d(round(size(frameUndistorted)*1.4));

    if currentFrame == 1
        % Only extract the coordinates of the buoy in the first frame.
        % GCF is the MATLAB key to the current figure.
        [xBuoy, yBuoy] = getInitialBuoyLocation(frame);
        framePrev = frameUndistorted;
    end
    
    if currentFrame >=2
        [frameUndistortedWarped, tform] = getCameraStabilisationTransform(frameUndistorted, framePrev, keyPointThreshold, worldMapping);
        
        % Memory for camera stabilisation    
        framePrev = frameUndistortedWarped;
        
        % ROI.
        frameRef = [xBuoy - 0.5*widthSearchArea yBuoy - 0.5*heightSearchArea];
        searchArea = frameUndistortedWarped(yBuoy - 0.5*heightSearchArea : yBuoy + 0.5*heightSearchArea,...
                                             xBuoy - 0.5*widthSearchArea : xBuoy + 0.5*widthSearchArea,...
                                             :);
        % Flow of searcharea
        flow = flowObj.estimateFlow(rgb2gray(searchArea));
                                         
        % Adaptive thresholding to retrieve a binarised image.
        Thres = adaptthresh(rgb2gray(searchArea), adaptThreshSensitivity);
        bin = imbinarize(rgb2gray(searchArea), Thres);
        
        % Blob analysis.
        [area, centroid, bbox, eccentricity, labeled] = blobInfo.step(bin);
        
        % Filter out blobs that have no flow or a flow that is too high.
        labeled2 = filterBasedOnFlow(labeled, flow.Magnitude, maxFlow);
        
        % Calculate probability that blob is the buoy.
        numberBlobs = size(eccentricity, 1);
        blobProb = zeros(numberBlobs,1);
        
        for b = 1:numberBlobs
            if sum(sum(labeled2 == b)) > 0
                blobProb(b,:) = (1/(eccentricity(b,:)+1))*h(round(centroid(b,1)), round(centroid(b,2)));     
            end
        end
        
        % Select the blob with the largest probability of being the buoy.
        [M,I] = max(blobProb(:));
        if (M >= minProb)
            xBuoy = round(centroid(I,1)) + frameRef(1);
            yBuoy = round(centroid(I,2)) + frameRef(2);
            autoLocationsBuoy(currentFrame, :) = [xBuoy yBuoy];
        else
            autoLocationsBuoy(currentFrame,:) = [-1 -1];
        end
       
        % Get the buoy coordinates in the unwarped, but undistorted frame.
        buoyOriginalImage = tform.transformPointsInverse([xBuoy yBuoy]);
        [buoyOriginalImage(1), buoyOriginalImage(2)] = worldToIntrinsic(worldMapping, buoyOriginalImage(1), buoyOriginalImage(2));
        
        % Determine the horizon.
        lines = horizonDetection(frame);
        
        % Interpolate the found horizon line over the entire width of the
        % image
        [rowsHorizon, colsHorizon] = interpolateHorizon(frame, lines);
        
        %% Calculate the distance to the buoy        
        % Use the distance to the horizon and the corresponding number of
        % pixels of the buoy to create an exponential base number. 
        % L2 is then the distance in pixels from the boat to the buoy along
        % the perpendicular vector to the horizon.
        % NOTE: swap commenting to enable alternate distance calculation.
        tempDistance = originalDistanceCalculation(frame, lines, rowsHorizon, buoyOriginalImage, realDistanceHorizon);
%         tempDistance = alternativeDistanceCalculation(lines,...
%                                                 cameraParams,...
%                                                 realDistanceHorizon,...
%                                                 buoyOriginalImage,...
%                                                 rowsHorizon);
        
        % Filter the distance calculations using a low-pass moving average
        % filter.
        distanceToBuoy(currentFrame) = tempDistance;
        
        % Actual low pass filter implementation.
        nrSample = 4;
        if (size(distanceToBuoy,2)>nrSample)
            distanceToBuoy(currentFrame) = 1/nrSample*(sum(distanceToBuoy(currentFrame-(nrSample-1):currentFrame)));
        end
%% Visualization.
        subplot(2,2,1)
        imshow(frameUndistortedWarped);
        hold on
        drawSearchGrid(xBuoy, yBuoy, widthSearchArea, heightSearchArea);
        rectangle('Position', [xBuoy-5, yBuoy-5, 10, 10], 'Curvature', [1 1], 'EdgeColor', 'g');
        hold off
        title('Original video feed with camera stabilisation');
        
        if writeOutputVideo
            % use the current warped frame and input search grid, circle
            % and buoy location.
            imageToWrite = insertMarker(frame, [buoyOriginalImage(1) buoyOriginalImage(2)], 's', 'Size', 50, 'Color', 'red');
            if autoLocationsBuoy(currentFrame,1)~=-1
                imageToWrite = insertMarker(imageToWrite, [buoyOriginalImage(1) buoyOriginalImage(2)], 'o', 'Size', 20, 'Color', 'green');
                imageToWrite = insertMarker(imageToWrite, [buoyOriginalImage(1) buoyOriginalImage(2)], 'x', 'Color', 'blue');
                imageToWrite = insertText(imageToWrite, [buoyOriginalImage(1)-50 buoyOriginalImage(2)-50], sprintf('Distance to target = %.0f m', distanceToBuoy(currentFrame)));
            end
            imageToWrite = insertText(imageToWrite,[50 50], sprintf('%d', currentFrame));
            outVideo.writeVideo(imageToWrite);
        end
        
        subplot(2,2,2)
        imshow(frame, [])
        hold on
        line([1 size(frame,2)], [rowsHorizon(1) rowsHorizon(size(frame,2))], 'Color', 'r', 'LineWidth', 1)
        drawSearchGrid(buoyOriginalImage(1), buoyOriginalImage(2), widthSearchArea, heightSearchArea);
        rectangle('Position', [buoyOriginalImage(1)-5, buoyOriginalImage(2)-5, 10, 10], 'Curvature', [1 1], 'EdgeColor', 'y');
        hold off
        title('Detected horizon and buoy in original frame');
        
        subplot(2,2,3)
        imshow(searchArea);
        hold on
        plot(widthSearchArea/2, heightSearchArea/2, 'r+');
        plot(flow);
        % Find quiver handle
        q = findobj(gca,'type','Quiver');
        % Change color to red
        q.Color = 'r';
        hold off;
        title('Zoomed searchgrid');
%         
        subplot(2,2,4) 
        plot(1:currentFrame, distanceToBuoy);
        title('Distance to buoy');
        xlabel('Frame number');
        ylabel('Distance (m)');
    end
    
    proc_time = toc
    drawnow limitrate
end

% Close the output video to actually make the data available outside
% matlab.
if writeOutputVideo
    close(outVideo);
end