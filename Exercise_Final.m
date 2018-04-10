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

%% Show the video
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
currentFrame = 0;

% Initialise initial buoy coordinates.
xBuoy = [];
yBuoy = [];

% Create optical flow object using Lucas-Kanade Derivative of Gaussian
flowObj = opticalFlowLKDoG( 'NoiseThreshold', 0.0039, 'NumFrames', 3,...
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

tform_translation = affine2d([1 0 0; 0 1 0; 50 50 1]);
                           
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
%             xBuoy = xBuoy+50; yBuoy=yBuoy+50;
            Buoy = tform_translation.transformPointsForward([xBuoy yBuoy]); % compensate for imtranslate further ahead
            xBuoy = Buoy(1); yBuoy = Buoy(2);
            frameUndistorted = imwarp(frameUndistorted, tform_translation);
        end
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
        frameUndistortedWarped = imwarp(frameUndistorted, tform, 'OutputView', imref2d(round(size(frameUndistorted)*1.4)));
        
        % Memory for camera stabilisation    
        framePrev = frameUndistortedWarped;   
        
%         frameUndistortedWarped = imtranslate(frameUndistortedWarped, [50 50]);
        frameUndistortedWarped = imwarp(frameUndistortedWarped, tform_translation, 'OutputView', imref2d(round(size(frameUndistorted)*1.4))); % translate the image to ensure four corners are always detected.
        
        % ROI.
        frameRef = [xBuoy - 0.5*widthSearchArea yBuoy - 0.5*heightSearchArea];
        searchArea = frameUndistortedWarped(yBuoy - 0.5*heightSearchArea : yBuoy + 0.5*heightSearchArea,...
                                             xBuoy - 0.5*widthSearchArea : xBuoy + 0.5*widthSearchArea,...
                                             :);
        % Flow of searcharea
        flow = flowObj.estimateFlow(rgb2gray(searchArea));
                                         
        % Morphological operations.
        Thres = adaptthresh(rgb2gray(searchArea), adaptThreshSensitivity);
        bin = imbinarize(rgb2gray(searchArea), Thres);
        
        % Blob analysis.
        [area, centroid, bbox, eccentricity, labeled] = blobInfo.step(bin);
        
        temp = logical(labeled).*flow.Magnitude;
        temp2 = temp>0;
        temp3 = temp<maxFlow;
        labeled2 = labeled.*uint8(temp3).*uint8(temp2);
%         bin2 = logical(temp4);
%         erod = imerode(bin2, strel('disk', 1));
%         dila = imdilate(erod, strel('disk', 1));
        
        % Calculate probability.
        numberBlobs = size(eccentricity, 1);
        blobProb = zeros(numberBlobs,1);
        
        for b = 1:numberBlobs
            if sum(sum(labeled2 == b)) > 0
                blobProb(b,:) = (1/(eccentricity(b,:)+1))*h(round(centroid(b,1)), round(centroid(b,2)));     
            end
        end
        
        [M,I] = max(blobProb(:));
        if (M >= minProb)
            xBuoy = round(centroid(I,1)) + frameRef(1);
            yBuoy = round(centroid(I,2)) + frameRef(2);                
        end
        temp = rgb2gray(frameUndistortedWarped) ~= 0;
        im_fx = ut_gauss(temp, 3, 1, 0);
        im_fy = ut_gauss(temp, 3, 0, 1);
        im_lap = im_fx + im_fy;
        cornerPoints = detectHarrisFeatures((im_lap));
        centers = subclust(cornerPoints.Location, 0.5);
        width_crop = inf; height_crop = inf;
        for i = 1:size(centers,1)
           for j = 1:size(centers,1)
              if i ~= j
                 temp_x = centers(i,1) - centers(j,1);
                 if  temp_x < width_crop && temp_x > 50 
                     width_crop = temp_x;
                     x_crop = centers(j,1);
                 end
                 temp_y = centers(i,2) - centers(j,2);
                 if  temp_y < height_crop && temp_y > 50
                     height_crop = temp_y;
                     y_crop = centers(j,2);
                 end
              end
           end
        end
        cropped = imcrop(frameUndistortedWarped, [x_crop+20,y_crop+20,width_crop-20,height_crop-20]);
        buoyCropped = [xBuoy - (x_crop+20), yBuoy - (y_crop+20)];
        edges = edge(rgb2gray(cropped), 'canny', [0.2 0.5]);
        [H, T, R] = hough(edges);
        P = houghpeaks(H,5);
        lines = houghlines(edges, T, R, P);
        colsHorizon = 1:size(cropped,2);
        rowsHorizon = (getfield(lines,{1}, 'point1',{2}) - getfield(lines,{1}, 'point2', {2})) / ... %dRow rows are second coordinate
                        (getfield(lines, {1}, 'point1', {1}) - getfield(lines, {1}, 'point2', {1})) * ... %dRow cols are first coordinate
                        colsHorizon;
        rowsHorizon = rowsHorizon + (getfield(lines,{1}, 'point1',{2}) - rowsHorizon(getfield(lines,{1}, 'point1',{1}))); % add the 'b' in ax+b to all values
        q_rad = deg2rad(-90 - getfield(lines, {1}, 'theta'));
        tform_horizonRot = affine2d([cos(q_rad) sin(q_rad) 0; -sin(q_rad) cos(q_rad) 0; 0 0 1]); %image doesn't rotate around center
        croppedHorizonCorrected = imwarp(cropped, tform_horizonRot, 'OutputView', imref2d(round(size(frameUndistorted)*1.4)));
        buoyCroppedRot = tform_horizonRot.transformPointsForward(buoyCropped); % compensation for rotation
        horizon = tform_horizonRot.transformPointsForward([colsHorizon' rowsHorizon']); % horizon is straigtend and equivalent to the new croppedHorizonCorrected
        nonLinearProfile = 5600^(1/(size(cropped,1)-mean(horizon(:,2)))); % 5600 distance to horizon in real life at 2.5 meter height.
        distanceToBuoy(currentFrame) = nonLinearProfile^(size(cropped,1)- buoyCroppedRot(2)); % needs filtering of NaNs, large differences, and fast transitions.
%% Visualization.
        subplot(2,2,1)
        imshow(frameUndistortedWarped);
        drawSearchGrid(xBuoy, yBuoy, widthSearchArea, heightSearchArea);
        rectangle('Position', [xBuoy-5, yBuoy-5, 10, 10], 'Curvature', [1 1], 'EdgeColor', 'g');
        title('Original video feed with camera stabilisation');
        
        if writeOutputVideo
            % use the current warped frame and input search grid, circle
            % and buoy location.
            imageToWrite = insertMarker(frameUndistortedWarped, Buoy, 's', 'Size', 100, 'Color', 'red');
            imageToWrite = insertMarker(imageToWrite, Buoy, 'o', 'Size', 20, 'Color', 'green');
            imageToWrite = insertMarker(imageToWrite, Buoy, 'x', 'Color', 'blue');
            imageToWrite = insertText(imageToWrite,[50 50], sprintf('%d', currentFrame));
            outVideo.writeVideo(imageToWrite);
        end
        
        subplot(2,2,2)
        imshow((croppedHorizonCorrected), [])
        hold on
%         line([1 size(cropped,2)], [rowsHorizon(1) rowsHorizon(size(cropped,2))], 'Color', 'r', 'LineWidth', 1)
%         max_len = 0;
%         for k = 1:length(lines)
%            xy = [lines(k).point1; lines(k).point2];
%            plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
% 
%            % Plot beginnings and ends of lines
%            plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%            plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
% 
%            % Determine the endpoints of the longest line segment
%            len = norm(lines(k).point1 - lines(k).point2);
%            if ( len > max_len)
%               max_len = len;
%               xy_long = xy;
%            end
%         end
        hold off
        title('Hough lines');
        
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
        
        subplot(2,2,4) 
        plot(1:currentFrame, distanceToBuoy);
%         imshow(label2rgb(labeled2));
%         hold on;
%         plot(widthSearchArea/2, heightSearchArea/2, 'r+');
%         hold off;
%         title('Filtered on threshold, area & flow');
    end
    
    proc_time = toc
    drawnow limitrate
end
if writeOutputVideo
    close(outVideo);
end