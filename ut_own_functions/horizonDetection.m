function lines = horizonDetection(frame)
% lines = HORIZONDETECTION(frame) detects the horizon by using the Houghtransform.
% FRAME is the image in which the horizon must be detected, and should be
% the unmodified video image.
%
% RETURNS the Hough lines for the horizon.
    edges = edge(rgb2gray(frame), 'canny', [0.2 0.5]);
    [H, T, R] = hough(edges);
    P = houghpeaks(H);
    lines = houghlines(edges, T, R, P);
end