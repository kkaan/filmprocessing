clear

% Default parameters
ROI_Pix_MLC = 1;
ROI_Pix_Norm = 5;
y_size = 768;
x_size = 1024;
PixelSize = 0.392;
SID = 100;
SDD = 180;
SaveResults = 1;

% Load Data
[Data_Mat_File,D] = uigetfile('*.mat','Select Data Mat File');
Data_Mat_File = [D Data_Mat_File(1:end-4)];
load(Data_Mat_File)

% Sort and read normalisation Data
NormImages = [];
MLCImages = [];
Num_Norm_Fields = numel(Norm_Data(:,1));
for i = 1:Num_Norm_Fields
    Energy_String = ['E_' char(Norm_Data{i,2})];
    Field_String = ['Field_' num2str(Norm_Data{i,1},'%d')];
    NormImages.(char(Energy_String)).(char(Field_String)).FileNums = Norm_Data(i,3:5);
    NormImages.(char(Energy_String)).(char(Field_String)).FieldSizemm = Norm_Data{i,1};
    Images = [];
    Header = [];
    for j = 1:3
        File = [Data_Dir 'SID0' num2str(Norm_Data{i,2+j}) '.dcm'];
        [Im,H,NFrames] = LoadEpidImage(File);
        Images{j} = Im*NFrames;
        Header{j} = H;
    end
    NormImages.(char(Energy_String)).(char(Field_String)).Images = Images;
    NormImages.(char(Energy_String)).(char(Field_String)).Header = Header;
end

% Sort and read mlc Data
Num_MLC_Fields = numel(MLC_Data(:,1));
for i = 1:Num_MLC_Fields
    Energy_String = ['E_' char(MLC_Data{i,2})];
    Field_String = ['Field_' num2str(MLC_Data{i,1},'%0.1f')];
    I = find(Field_String == '.');
    Field_String(I) = 'p';
    MLCImages.(char(Energy_String)).(char(Field_String)).FileNums = MLC_Data(i,3:5);
    MLCImages.(char(Energy_String)).(char(Field_String)).FieldSizemm = MLC_Data{i,1};
    Images = [];
    Header = [];
    for j = 1:3
        File = [Data_Dir 'SID0' num2str(MLC_Data{i,2+j}) '.dcm']; % '%05d'
        [Im,H,NFrames] = LoadEpidImage(File);
        Images{j} = Im*NFrames;
        Header{j} = H;
    end
    MLCImages.(char(Energy_String)).(char(Field_String)).Images = Images;
    MLCImages.(char(Energy_String)).(char(Field_String)).Header = Header;
end

% Analyse Now energy by energy
Norm_Energies = fieldnames(NormImages);
Num_Energies = numel(Norm_Energies);

for E = 1:Num_Energies
    Energy_String = Norm_Energies{E};
    % Get CAX for Norm Fields
    Norm_Mean = [];
    Jaw_OR_Mean =[];
    Jaw_OR_Std = [];
    Jaw_FieldSize = [];
    Norm_Field_Names = fieldnames(NormImages.(char(Energy_String)));
    Num = numel(Norm_Field_Names);
    for i = 1:Num
        Images = NormImages.(char(Energy_String)).(char(Norm_Field_Names{i})).Images;
        for j = 1:3
            Im = Images{j};
            [X_CAX, Y_CAX] = Get_CAX_Cones(Im);
            [Mean_Pix(j),Std_Pix(j)] = Get_CAXPix_Cones(Im,X_CAX,Y_CAX,ROI_Pix_Norm);
        end
        Norm_Mean.(char(Norm_Field_Names{i})) = mean(Mean_Pix);
        Jaw_OR_Mean(i) = mean(Mean_Pix)./(Norm_Mean.Field_100);
        Jaw_OR_Std(i) = std(Mean_Pix)./(Norm_Mean.Field_100);
        Jaw_FieldSize(i) = NormImages.(char(Energy_String)).(char(Norm_Field_Names{i})).FieldSizemm;
    end
    % Get CAX for MLC Fields
    MLC_OR_Mean =[];
    MLC_OR_Std = [];
    MLC_FieldSize = [];
    MLC_X_Mean = [];
    MLC_Y_Mean = [];
    MLC_X_Std = [];
    MLC_Y_Std = [];
    MLC_Field_Names = fieldnames(MLCImages.(char(Energy_String)));
    Num = numel(MLC_Field_Names);
    for i = 1:Num
        Images = MLCImages.(char(Energy_String)).(char(MLC_Field_Names{i})).Images;
        X_Profile = zeros(3,1024);
        Y_Profile = zeros(3,768);
        for j = 1:3
            Im = Images{j};
            [X_CAX(j), Y_CAX(j)] = Get_CAX_Cones(Im);
            [Mean_Pix(j),Std_Pix(j)] = Get_CAXPix_Cones(Im,X_CAX(j),Y_CAX(j),ROI_Pix_MLC);
            [X_mm, Y_mm, X_Profile(j,:), Y_Profile(j,:)] = Get_Profiles_Cones(Im,X_CAX(j),Y_CAX(j),SDD,SID,PixelSize);
        end
        
        Data.(char(Energy_String)).(char(MLC_Field_Names{i})).measurementVals = Mean_Pix;
        Data.(char(Energy_String)).(char(MLC_Field_Names{i})).CVMeasurements = std(Mean_Pix)/mean(Mean_Pix);
        Data.(char(Energy_String)).(char(MLC_Field_Names{i})).meanMeasurements = mean(Mean_Pix);
  
        MLC_OR_Mean(i) = mean(Mean_Pix)./(Norm_Mean.Field_100);
        MLC_OR_Std(i) = std(Mean_Pix)./(Norm_Mean.Field_100);
        MLC_FieldSize(i) = MLCImages.(char(Energy_String)).(char(MLC_Field_Names{i})).FieldSizemm;
        MLC_X_Mean(i) = (mean(X_CAX(j)) - x_size/2)*PixelSize*SID/SDD;
        MLC_X_Std(i) = std(X_CAX(j))*PixelSize*SID/SDD;
        MLC_Y_Mean(i) = (mean(Y_CAX(j)) - y_size/2)*PixelSize*SID/SDD;
        MLC_Y_Std(i) = std(Y_CAX(j))*PixelSize*SID/SDD;
        MLC_X_Coord{i} = X_mm';
        MLC_Y_Coord{i} = Y_mm';
        MLC_X_Profile{i} = mean(X_Profile)/max(mean(X_Profile))*MLC_OR_Mean(i);
        MLC_Y_Profile{i} = mean(Y_Profile)/max(mean(Y_Profile))*MLC_OR_Mean(i);
    end
    if strcmp(Linac,'G1')
        C = 'k';
    else
        C = 'r';
    end
%     C = 'b';
    % Plot Output Ratios
    h_OR = figure;
    errorbar(MLC_FieldSize,MLC_OR_Mean,MLC_OR_Std,C,'LineWidth',2)
    title(['Output Ratio - ' Energy_String(3:end)])
    xlabel('Field Size (mm)')
    ylabel('Output Ratio')
    grid on
    % Plot X and Y Centres
    h_CAX = figure;
    subplot(2,1,1)
    plot(MLC_FieldSize,MLC_X_Mean - mean(MLC_X_Mean),char(['o' C]),'MarkerSize',5,'LineWidth',2)
    title(['Concentricty X Direction - ' Energy_String(3:end)])
    xlabel('Field Size (mm)')
    ylabel('Field Offset (mm)')
    axis([4 21 -0.3 +0.3])
    grid on
    subplot(2,1,2)
    plot(MLC_FieldSize,MLC_Y_Mean - mean(MLC_Y_Mean),char([C 'o']),'MarkerSize',5,'LineWidth',2)
    title(['Concentricty Y Direction - ' Energy_String(3:end)])
    xlabel('Cone Diameter (mm)')
    ylabel('Cone Offset (mm)')
    axis([4 21 -0.3 +0.3])
    grid on    
    % Plot Profiles
    h_Prof = figure;
    for i = 1:Num
        subplot(2,1,1)
        plot(MLC_X_Coord{i},MLC_X_Profile{i},C,'LineWidth',2)
        title(['X Profiles - ' Energy_String(3:end)])
        xlabel('Distance from Centre of Field (mm)')
        ylabel('Normalised Output')
        hold on
        grid on
        axis([-20 20 0 1])
        subplot(2,1,2)
        plot(MLC_Y_Coord{i},MLC_Y_Profile{i},C,'LineWidth',2)
        title(['Y Profiles - ' Energy_String(3:end)])
        xlabel('Distance from Centre of Field (mm)')
        ylabel('Normalised Output')
        grid on
        hold on
        axis([-20 20 0 1])
    end
    if SaveResults
        savefig(h_OR,[Data_Mat_File '_' Energy_String(3:end) '_MLC_Output Ratio.fig'])
        savefig(h_CAX,[Data_Mat_File '_' Energy_String(3:end) '_MLC_Concentricity.fig'])
%        savefig(h_Prof,[Data_Mat_File '_' Energy_String(3:end) '_MLC_Profiles.fig'])
        SaveFile = [Data_Mat_File '_' Energy_String(3:end) '_MLC_Results'];
        save(SaveFile,'MLC_FieldSize','MLC_OR_Mean','MLC_OR_Std','MLC_X_Mean','MLC_Y_Mean','MLC_X_Profile','MLC_Y_Profile');
    end
end
