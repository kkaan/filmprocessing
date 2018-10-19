function [X_mm, Y_mm, X_Profile, Y_Profile] = filmprofiles(Image, X_CAX, Y_CAX, pixelsize)
    
% Set Debug Value
DebugMode = 0;

% How many profiles to include on either side of CAX profile?
Num = 1;

% Get the image dimensions
[y_size, x_size] = size(Image);

% Get X and Y Pix
X_Pix = round(X_CAX);
Y_Pix = round(Y_CAX);

% Get profiles from image
Y_Profiles_Raw = Image(:,X_Pix-Num:X_Pix+Num);
X_Profiles_Raw = Image(Y_Pix-Num:Y_Pix+Num,:);
Y_Profiles_Raw = Y_Profiles_Raw';

% Average profiles
X_Profiles_Mean = mean(X_Profiles_Raw,1);
Y_Profiles_Mean = mean(Y_Profiles_Raw,1);

% Do some smoothing
X_Profile = X_Profiles_Mean;
Y_Profile = Y_Profiles_Mean;

% Get the X and Y coordinate vectors in mm
X_mm = ((1:x_size) - X_CAX)*pixelsize;
Y_mm = ((1:y_size) - Y_CAX)*pixelsize;

% PLOTS TO CHECK
if DebugMode
    figure
    plot(X_mm,X_Profile,'b','Linewidth',2)
    hold on
    plot(Y_mm,Y_Profile,'r','Linewidth',2)
    grid on
    xlabel('Distance at Isocentre (mm)')
    ylabel('Raw Pixel Values (arb. units)')
    legend('X Profile','Y Profile')
    set(gca,'XLim',[-20 20])
end