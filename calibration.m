function [NetOD, p] = calibration(calibrationFiles, path, Doses, polyorder)

% Takes in a list of filenames of images to create a calibration curve
% 
if (nargin < 1)
    [calibrationFiles, path, ~] = uigetfile('*.tiff', 'Select calibration images', 'MultiSelect', 'On');
    Doses = [0, 75, 150, 225, 300, 375, 450, 525, 600];
    polyorder = 2;
end



% Check is number of calibrations match the number of doses provided.

% For each calibration image obtain the average pixel value of the central
% quadrant. Save to array
NetOD = zeros(1, numel(calibrationFiles));


for n = 1:numel(calibrationFiles)
    netODImageFile = fullfile(path,calibrationFiles{n,1});
    netODIm = imread(netODImageFile);
    
    % Crop to the central quarter of the image
    
    horizmargin     = int16(size(netODIm,1)/3);
    vertmargin      = int16(size(netODIm,2)/3);
    
    
    netODCropped    = netODIm(horizmargin:end-horizmargin ,...
                        vertmargin:end-vertmargin);
    NetOD(1,n)      = mean(mean(netODCropped));
 
end

% fit a third degree polynomial to the curve.
p = polyfit(NetOD, Doses, polyorder); 



% plot(NetOD, Doses, 'o')
% x1 = [0:0.01:max(NetOD)+0.2]
% y1 = polyval(p, x1)
% 
% figure
% plot(x,y,'o')
% hold on
% plot(x1,y1)
% hold off



