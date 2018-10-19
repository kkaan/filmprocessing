function FOFcone(Norm_Data, Cone_Data, LinacName, Meas_Date, Data_Dir)
% Produce output factors and profiles from calibrated film dose files


% Default parameters
ROI_Pix_Cone = 4;
ROI_Pix_Norm = 15;
dpi=75;
PixelSize = 25.4/dpi;
SID = 100;
SDD = 180;
SaveResults = 1;

prefix = '';
cropmargin = 20;

Norm_MU = 500;
Field_MU   = 700;


% Sort and read normalisation Data
NormImages = [];
ConeImages = [];
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

% Sort and read cone Data
Num_Cone_Fields = numel(Cone_Data(:,1));
for i = 1:Num_Cone_Fields
    Energy_String = ['E_' char(Cone_Data{i,2})];
    Field_String = ['Field_' num2str(Cone_Data{i,1},'%0.1f')];
    I = find(Field_String == '.');
    Field_String(I) = 'p';
    ConeImages.(char(Energy_String)).(char(Field_String)).FileNums = Cone_Data(i,4:6);
    ConeImages.(char(Energy_String)).(char(Field_String)).FieldSizemm = Cone_Data{i,1};
    im = [];
    Header = [];
    nim = sum(~isnan([Cone_Data{i, 4:6}]));
    for j = 1:nim
        file        = fullfile(Data_Dir, [num2str(Cone_Data{i,3+j}) '.tiff']);
        im{j}       = imread(file{1});
    end
    ConeImages.(char(Energy_String)).(char(Field_String)).Images = im;
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
    Mean_Pix = zeros(1,Num); %reset between each cone
    Std_Pix  = zeros(1,Num); %reset between each cone
    
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
    % Get CAX for Cone Fields
    Cone_OR_Mean =[];
    Cone_OR_Std = [];
    Cone_FieldSize = [];
    Cone_X_Mean = [];
    Cone_Y_Mean = [];
    Cone_X_Std = [];
    Cone_Y_Std = [];
    Prof_all = [];
    cross_Prof = [];
    cross_fs = [];
    crossplane_av = [];
    x_pr = [];
    
    Cone_Field_Names = fieldnames(ConeImages.(char(Energy_String)));
    Num = numel(Cone_Field_Names);
    
        % set variables for resampling/cropping of film cross-profile data  
        x = [];
        y = [];
        diffx  = [];
        diffx2  = [];
        x_mm_new = [];
        x_pr_new = [];
        new_start = -15;
        new_end = 15;
        fs = 5;        
        startcrop = abs(-17-new_start)*fs+1;
        endcrop = (17-new_end)*fs;
        film_len = (new_end - new_start)*fs+1;
        xx_mm = [];
        xx_pr = [];
    
    for i = 1:Num
        im = ConeImages.(char(Energy_String)).(char(Cone_Field_Names{i})).Images;
        [~,nim]                 = size(im);
        Im                      = im{1};
        imcr                    = Im(cropmargin:end-cropmargin, cropmargin:end-cropmargin);
        
        [Y_range,X_range]       = size(imcr);
        X_Profile               = zeros(nim, X_range);
        X_mm                    = zeros(nim, X_range);
        Y_Profile               = zeros(nim, Y_range);
        Y_mm                    = zeros(nim, Y_range);
        
        Mean_Pix = zeros(1,nim); %reset between each cone
        Std_Pix  = zeros(1,nim); %reset between each cone
        
        Prof_all{i} = zeros(1,nim); %reset between each cone
        cross_Prof{i} = zeros(1,nim); %reset between each cone
        cross_fs{i} = zeros(1,nim);
        crossplane_av{i} = zeros(film_len,1); %reset between each cone
                
        crossplane = figure(1);
        title('Crossplane profile')
        xlabel('Off-axis distance (mm)') % x-axis label
        ylabel('Dose (cGy)') % y-axis label
        
            
        inplane = figure(2);
        title('In-plane profile')
        xlabel('Off-axis distance (mm)') % x-axis label
        ylabel('Dose (cGy)') % y-axis label
   
        
        
        for j = 1:nim                    
            Im                          = im{j};
            % crop margins off to get rid of bad bits of film
            Im                          = Im(cropmargin:end-cropmargin, cropmargin:end-cropmargin);
            [X_CAX(j), Y_CAX(j)]        = Get_CAX_Cones(Im);
            [Mean_Pix(j),Std_Pix(j)]    = Get_CAXPix_Cones(Im,X_CAX(j),Y_CAX(j),ROI_Pix_Cone);
            [X_mm(j,:), Y_mm(j,:), X_Profile(j,:), Y_Profile(j,:)] ...
                                        = filmprofiles(Im,X_CAX(j),Y_CAX(j),PixelSize);
            
            %% resample cros-profile film data so that all profiles have common distance axes:
            a = size(X_mm(j,:));
            xlen = a(2);
                                    
            %fill film data to span -17 mm to 17 mm
            diffx{j} = -17-X_mm(j,1);       % extend profile data to start at -17cm
            diffx2{j,1} = 17-X_mm(j,xlen);  % extend profile data to end at +17cm

            x_mm_new{j,1}   = [diffx{j}+X_mm(j,1)  X_mm(j,:)        diffx2{j} + X_mm(j,xlen)];
            x_pr_new{j}   =   [X_Profile(j,1)      X_Profile(j,:)   X_Profile(j,xlen)];

            %resample new film data, and crop down to -15 mm to 15 mm
            [X_pr1{j}, x1{j}] = resample(x_pr_new{j},x_mm_new{j},fs,3,1);
            xx_pr{j} = X_pr1{j}(startcrop:end-endcrop);
            xx_mm{j} = x1{j}(startcrop:end-endcrop);
 
            cross_fs{i} = xx_mm{j}'; %cross profile field size
            x_pr{i,j}   = xx_pr{j}'/max(xx_pr{j});% normalised cross profile data for each measurement
            
            %%
            figure(1)
            hold on
            % plot(X_mm(j,:), 100*X_Profile(j,:)/max(X_Profile(j,:)));
            plot(xx_mm{j},100*(xx_pr{j}/max(xx_pr{j})))
            
            figure(2)
            hold on
            plot(Y_mm(j,:), Y_Profile(j,:));
                    
        end
   
        hold off
        
        Cone_OR_Mean(i) = (mean(Mean_Pix)/Field_MU)./((Norm_Mean.Field_100)/Norm_MU);
        Cone_OR_Std(i) = std(Mean_Pix)./(Norm_Mean.Field_100);
        Cone_FieldSize(i) = ConeImages.(char(Energy_String)).(char(Cone_Field_Names{i})).FieldSizemm;
        
    
        

    end
    if strcmp(LinacName{1},'Troy')
        C = 'r';
    else
        C = 'b';
    end
%     C = 'b';
    % Plot Output Ratios
    h_OR = figure;
    errorbar(Cone_FieldSize,Cone_OR_Mean,Cone_OR_Std,C,'LineWidth',2)
    title(['Output Ratio - ' Energy_String(3:end)])
    xlabel('Cone Diameter (mm)')
    ylabel('Output Ratio')
    grid on

    X_PR{1} = [[x_pr{1,1}] [x_pr{1,2}] [x_pr{1,3}]];
    X_PR{2} = [[x_pr{2,1}] [x_pr{2,2}] [x_pr{2,3}]];
    X_PR{3} = [[x_pr{3,1}] [x_pr{3,2}] [x_pr{3,3}]];
    X_PR{4} = [[x_pr{4,1}] [x_pr{4,2}] [x_pr{4,3}]];
    X_PR{5} = [[x_pr{5,1}] [x_pr{5,2}]];
    X_PR{6} = [[x_pr{6,1}] [x_pr{6,2}]];
    X_PR{7} = [[x_pr{7,1}] [x_pr{7,2}]];
    for i = 1:Num   
        X_PR_mean{i} = mean(X_PR{i},2);
    end

    if SaveResults
        
        crsplfile = [LinacName{1} Energy_String(3:end) '_crossplane.fig'];
        crsplfilefull = fullfile(Data_Dir{1}, crsplfile);
        savefig(crossplane, crsplfilefull)
        
        inplfile = [LinacName{1} Energy_String(3:end) '_inplane.fig'];
        inplfilefull = fullfile(Data_Dir{1}, inplfile);
        savefig(inplane,inplfilefull);
        
        orfile = [LinacName{1} '_' Energy_String(3:end) '_output_ratios.fig'];
        orfilefull = fullfile(Data_Dir{1}, orfile);
        savefig(h_OR,orfilefull)
        
        SaveFile = [LinacName{1} '_' Energy_String(3:end) '_Results'];
        SaveFileFull = fullfile(Data_Dir{1}, SaveFile);
        save(SaveFileFull,'Cone_FieldSize','Cone_OR_Mean','Cone_OR_Std','x_pr','cross_fs','X_PR_mean');
    end
end

%end of function 
end


%   X_std = std(X,0, 2);
%   fwhm(xx{1},X_mean)
%   fwhm(xx{1},X_mean+2*X_std)
%   fwhm(xx{1},X_mean-2*X_std)

