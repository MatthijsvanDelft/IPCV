function drawSearchGrid(xPos, yPos, width, height)
%         % Draw the search grid in the image
        rectangle( 'Position',[xPos-0.5*width,...
                   yPos-0.5*height, width,...
                   height], 'EdgeColor', 'r', 'LineWidth', 1,...
                   'LineStyle','-', 'Curvature', 0.2)
end