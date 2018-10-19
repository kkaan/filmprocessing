%function FOFmlc(Norm_Data, MLC_Data, LinacName, Meas_Date, Data_Dir)
% Produce output factors and profiles from calibrated film dose files
clear


% Default parameters
ROI_Pix_MLC = 4;
ROI_Pix_Norm = 15;
dpi=75;
PixelSize = 25.4/dpi;
SID = 100;
SDD = 180;
SaveResults = 1;

prefix = '';
cropmargin = 20;
cropmargin_1 = 10;

Norm_MU = 500;
Field_MU   = 600;

% Load Data
 [excelfile, path, ~] = uigetfile({ '*.*', 'All Files'},'Select excel file describing film measurements');
 excelfile = fullfile(path, excelfile);
 [Cone_Data, MLC_Data, LinacName, Meas_Date, Norm_Data, Data_Dir] = filmexcelextract(excelfile);

% Sort and read normalisation Data
NormImages = [];
MLCImages = [];
Num_Norm_Fields = numel(Norm_Data(:,1));
for i = 1:Num_Norm_Fields
    Energy_String = ['E_' char(Norm_Data{i,2})];
    Field_String = ['Field_' num2str(Norm_Data{i,1},'%d')];
    NormImages.(char(Energy_String)).(char(Field_String)).FileNums = Norm_Data(i,4:6);
    NormImages.(char(Energy_String)).(char(Field_String)).FieldSizemm = Norm_Data{i,1};
    NormImages.(char(Energy_String)).(char(Field_String)).MU = Norm_Data{i,3};
    im = [];
    Header = [];
    nim = sum(~isnan([Norm_Data{i, 4:6}]));
    for j = 1:nim
        if isnumeric(Norm_Data{i,3+j})
            file        = fullfile(Data_Dir, [num2str(Norm_Data{i,3+j}) '.tiff']);
        else
            file        = fullfile(Data_Dir,[Norm_Data{i,3+j} '.tiff']);
        end
        im{j}       = imread(file{1});
    end
    NormImages.(char(Energy_String)).(char(Field_String)).Images = im;
end

% Sort and read mlc Data
Num_MLC_Fields = numel(MLC_Data(:,1));
for i = 1:Num_MLC_Fields
    Energy_String = ['E_' char(MLC_Data{i,2})];
    Field_String = ['Field_' num2str(MLC_Data{i,1},'%0.1f')];
    I = find(Field_String == '.');
    Field_String(I) = 'p';
    MLCImages.(char(Energy_String)).(char(Field_String)).FileNums = MLC_Data(i,4:6);
    MLCImages.(char(Energy_String)).(char(Field_String)).FieldSizemm = MLC_Data{i,1};
    im = [];
    Header = [];
    nim = sum(~isnan([MLC_Data{i, 4:6}]));
    for j = 1:nim
        file        = fullfile(Data_Dir, [num2str(MLC_Data{i,3+j}) '.tiff']);
        im{j}       = imread(file{1});
    end
    MLCImages.(char(Energy_String)).(char(Field_String)).Images = im;
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
    Mean_Pix = zeros(1,Num); %reset between each field size
    Std_Pix  = zeros(1,Num); %reset between each field size
    
    for i = 1:Num
        im = NormImages.(char(Energy_String)).(char(Norm_Field_Names{i})).Images;
        [~,nim] = size(im);
        for j = 1:nim
            Im                  = im{j};
            Im                  = Im(cropmargin:end-cropmargin, cropmargin:end-cropmargin);
            [y_size, x_size]    = size(Im);
            Y_CAX               = round(y_size/2);
            X_CAX               = round(x_size/2); 
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
        [~,nim]                 = size(Images);
        Im                      = Images{1};
        imcr                    = Im(cropmargin:end-cropmargin, cropmargin:end-cropmargin);
        
        [Y_range,X_range]       = size(imcr);
        X_Profile               = zeros(nim, X_range);
        X_mm                    = zeros(nim, X_range);
        Y_Profile               = zeros(nim, Y_range);
        Y_mm                    = zeros(nim, Y_range);
        
        for j = 1:3
            Im = Images{j};
            Im = Im(cropmargin_1:end-cropmargin_1, cropmargin_1:end-cropmargin_1);
            [y_size, x_size] = size(Im);
     
            X_Profile = zeros(3,x_size);
            Y_Profile = zeros(3,y_size);
            
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
    
 %   if strcmp(LinacName{1},'Troy')
        C = 'r';
%    else
        C = 'b';
%    end
%     C = 'b';
    % Plot Output Ratios
    h_OR = figure;
    errorbar(MLC_FieldSize,MLC_OR_Mean,MLC_OR_Std,C,'LineWidth',2)
    title(['Output Ratio - ' Energy_String(3:end)])
    xlabel('MLC Diameter (mm)')
    ylabel('Output Ratio')
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
        
%        crsplfile = [LinacName{1} Energy_String(3:end) '_crossplane.fig'];
%        crsplfilefull = fullfile(Data_Dir{1}, crsplfile);
%        savefig(crossplane, crsplfilefull)
        
%        inplfile = [LinacName{1} Energy_String(3:end) '_inplane.fig'];
 %       inplfilefull = fullfile(Data_Dir{1}, inplfile);
%        savefig(inplane,inplfilefull);
        
%        orfile = [LinacName{1} '_' Energy_String(3:end) '_output_ratios.fig'];
%        orfilefull = fullfile(Data_Dir{1}, orfile);
%        savefig(h_OR,orfilefull)
        
%        SaveFile = [LinacName{1} '_' Energy_String(3:end) '_Results'];
%        SaveFileFull = fullfile(Data_Dir{1}, SaveFile);
%        save(SaveFileFull,'MLC_FieldSize','MLC_OR_Mean','MLC_OR_Std');

        savefig(h_Prof,[excelfile '_' Energy_String(3:end) '_MLC_Profiles.fig'])
        SaveFile = [excelfile '_' Energy_String(3:end) '_MLC_Results.mat'];
        save(SaveFile,'MLC_FieldSize','MLC_OR_Mean','MLC_OR_Std','MLC_X_Mean','MLC_Y_Mean','MLC_X_Profile','MLC_Y_Profile');
     end
end

%end of function 
%end

