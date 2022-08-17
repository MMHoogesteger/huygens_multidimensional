function [checkers] = findCheckers(camIm)
%findCheckers(camIm) Find Checkerboards in the image
% Finds all checkerboards in the image, check if boars are valid and of the
% correct size. To accomplish this, once a checkerbaord is found, all white
% checkers are made black and the checkerboard algorithm is run again.
%
% See also: annotateFromFindImage, findMarkers


% Convert the image to the HSV color space.
imHSV = rgb2hsv(camIm);

% Get the saturation channel.
value = imHSV(:, :, 3);

% Threshold the image
t = graythresh(value);
imBW = (value > t);

%% Checkerboards
hasCheckerboard = true;
c = 0;
checkers = struct;
checkers.n =0;
while(hasCheckerboard)
    [imagePoints, boardSize] = detectCheckerboardPoints(imBW);
    validBoard = ~isempty(imagePoints);
    validBoard = validBoard && all(abs(max(imagePoints)-min(imagePoints))>5);
    if(validBoard)
        wPoint = (imagePoints(2,:) + imagePoints(3,:)...
                + imagePoints(boardSize(1),:) - imagePoints(1,:))./2;
        bPoint = (imagePoints(1,:) + imagePoints(2,:)...
                + imagePoints(boardSize(1),:) - imagePoints(1,:))./2;
        wPointR = round(wPoint);
        bPointR = round(bPoint);
        if(~imBW(wPointR(2),wPointR(1)))
            if(~imBW(bPointR(2),bPointR(1)))
                % Both white and black point are black, stop
                hasCheckerboard = false;
            else
                wPoint = bPoint;
            end
        end

        tolerance = 1;
        whiteChecks = grayconnected(uint8(imBW*255), round(wPoint(2)), round(wPoint(1)), tolerance);

        imBW(whiteChecks) = 0;
        if(xor(mod(boardSize(1),2),mod(boardSize(2),2)))
            c = c+1;

            checkers(c).imagePoints = imagePoints;
            checkers(c).boardSize = boardSize;
            checkers(c).wPoint = wPoint;
            checkers(c).whiteChecks = whiteChecks;

        end
    else
        hasCheckerboard = false;
    end

end

for idx = 1:c
    checkers(idx).n = c;
end

