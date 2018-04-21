function [frameWarped, tform] = getCameraStabilisationTransform(frame, framePrev, keyPointThreshold, worldMapping)
% [frameWarped, tform] = GETCAMERASTABILISATIONTRANSFORM(frame, framePrev, keyPointThreshold, worldMapping)
% provides functionality to perform camera stabilisation given two input images, FRAME and FRAMEPREV. FRAME is the
% current image that needs to be stabilised, and FRAMEPREV is the previous
% stabilised image. 
%
% This algorithm uses FAST features to detect keypoints in the images and
% requires KEYPOINTTHRESHOLD which is the minimum intensity difference
% betwen corner and surrounding region, specified as the comma-seperated.
% The scalar value must be in the range (0,1). WORLDMAPPING is the of the
% type AFFINE2D and indicates the world limits into which the stablised
% image is to be warped.
%
% RETURNS the stabilised image, FRAMEWARPED and the geometrical transform
% used to stabilise the image, TFORM.
    % Keypoint detection.
    pointsCur = detectFASTFeatures(rgb2gray(frame), 'MinContrast', keyPointThreshold);
    pointsPrev = detectFASTFeatures(rgb2gray(framePrev), 'MinContrast', keyPointThreshold);

    % Feature extraction.
    [featuresCur, pointsCur] = extractFeatures(rgb2gray(frame), pointsCur);
    [featuresPrev, pointsPrev] = extractFeatures(rgb2gray(framePrev), pointsPrev);

    % Keypoint matching.
    indexPairs = matchFeatures(featuresPrev, featuresCur);
    pointsPrev = pointsPrev(indexPairs(:,1), :);
    pointsCur = pointsCur(indexPairs(:,2), :);

%         showMatchedFeatures(framePrev, frame, pointsPrev, pointsCur);

    % Stabilization transformation.
    [tform, ~, ~] = estimateGeometricTransform(pointsCur, pointsPrev, 'similarity');
    [frameWarped,~] = imwarp(frame, tform, 'OutputView', worldMapping);
end