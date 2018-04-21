function [rows, cols] = interpolateHorizon(frame, lines)
% [rows, cols] = INTERPOLATEHORIZON(frame, lines) interpolates the given
% LINES over the given FRAME. In this project, the line will be the
% horizon.
%
% RETURNS the column and row indices of the horizon in the entire width of
% the frame.
    cols = 1:size(frame,2);
    rowsHorizon = (getfield(lines,{1}, 'point1',{2}) - getfield(lines,{1}, 'point2', {2})) / ... %dRow rows are second coordinate
                    (getfield(lines, {1}, 'point1', {1}) - getfield(lines, {1}, 'point2', {1})) * ... %dRow cols are first coordinate
                    cols;
    rows = rowsHorizon + (getfield(lines,{1}, 'point1',{2}) - rowsHorizon(getfield(lines,{1}, 'point1',{1}))); % add the 'b' in ax+b to all values
end