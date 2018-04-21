function dist = originalDistanceCalculation(frame, lines, rowsHorizon, buoyOriginalImage, realDistanceHorizon)
    % Use the found information from the hough transform, and use
    % trigonometry to determine the perpendicular distance from the
    % boat to the horizon.
    phi = deg2rad(90 - abs(getfield(lines, {1}, 'theta')));
    temp_r1 =  buoyOriginalImage(2) - rowsHorizon(round(buoyOriginalImage(2)));
    L1 = (temp_r1)*cos(phi);
    L2 = (size(frame,1) - buoyOriginalImage(2)) / cos(phi);
    d_horizon = L1 + L2;
    dist = realDistanceHorizon^(L2/d_horizon);
end