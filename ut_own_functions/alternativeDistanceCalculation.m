function dist = alternativeDistanceCalculation(lines, cameraParams, realDistanceHorizon, buoyOriginalImage, rowsHorizon)
% dist = ALTERNATIVEDISTANCECALCULATION(lines, cameraParams,
% realDistanceHorizon, buoyOritinalImage, rowsHorizon) calcualates the
% distance to the BUOYORIGINALIMAGE using the angles found by the Hough
% transform, given by LINES, the camera parameters, CAMERAPARAMS, the
% distance to the actual horizon, REALDISTANCEHORIZON, and all the row
% indices of the horizon, ROWSHORIZON.
%
% RETURNS the distance to the buoy according to the alternative method.
    % Create the rotation matrix to realign the horizon.
    q_rad = deg2rad(-90 - getfield(lines, {1}, 'theta'));

    % Calculate the rotation to align principle point with horizon
    YTranslation = cameraParams.PrincipalPoint(2) - rowsHorizon(round(cameraParams.PrincipalPoint(1)));
    theta = tan(abs(YTranslation)/cameraParams.FocalLength(2));

    % Create ZXY Euler rotation matrix to correct the horizon (making
    % it horizontal) and aligning the principle axis with the horizon.
    eul = [theta 0 q_rad];
    rotmZYX = eul2rotm(eul);

    % Calculate the distance using camera parameters and euler matrix.
%         cRw = eye(3);
    cRw = rotmZYX;
    ctw = [0;0;realDistanceHorizon];
    cMat = [cRw, ctw; zeros(1,3), 1];
    K = cameraParams.IntrinsicMatrix';
    M = K * [cRw ctw];
%         M = K * [eye(3) zeros(3,1)];
%         M = M * cMat;
    p = [buoyOriginalImage(1); buoyOriginalImage(2); 1];
    wX = M\p;
    
    dist = wX(3);
end