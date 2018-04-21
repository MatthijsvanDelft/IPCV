function [x, y] = getInitialBuoyLocation(frame)
% GETINITIALBUOYLOCATION uses getpts to let the user select the inital position of the buoy.
% Requires an active figure window to show frame and to use getpts.
% 
% RETURNS x and y, in column and row indices respectively. Unit is pixels.
    xBuoy = [];
    yBuoy = [];
    
    imshow(frame);
    while (size(xBuoy, 1) > 1 || size(yBuoy, 1) > 1) || (isempty(xBuoy) || isempty(yBuoy))
        uiwait(msgbox({'Please pick one point.';...
                       'A shift-, right-, or double-click adds a final point and ends the selection.';...
                       'Pressing Return or Enter ends the selection without adding a final point. ';...
                       'Pressing Backspace or Delete removes the previously selected point.'}, 'Attention!', 'help'));
        [xBuoy, yBuoy] = getpts(gcf);
    end
    x = xBuoy;
    y = yBuoy;
end