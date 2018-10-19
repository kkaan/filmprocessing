function file = saveCal2Txt(p, doses, doseFiles, dosePath)
%SAVECAL2TXT saves the calibration coefficients and details into 
%text file

%   Detailed explanation goes here
    coeffstrr = sprintf('%.5f ' , p);
    dosestrr = sprintf('%.2f ', doses);
    imagefilestrr = sprintf('%s, ', doseFiles{:});
    
    %create the dose file folder
    idcs            = strfind(dosePath,'\');
    parentdir       = dosePath(1:idcs(end-1));
    
%     prompt = {'Energy:','Date:', 'Max Dose:'};
%     title = 'Calibration Param';
%     dims = [1 35];
%     definput = {'20','hsv'};
%     answer = inputdlg(prompt,title,dims,definput)
%     
    % YYMMDD_EXX_MAXXXXX_Channel_Cal.txt
    
    file = fullfile(parentdir, 'caliFileName.txt');
    fid = fopen(file, 'w');
    fprintf(fid, 'Polynomial Coefficients: ');
    fprintf(fid, coeffstrr);
    fprintf(fid, '\r\n\r\n');
    fprintf(fid, 'Doses: ');
    fprintf(fid, dosestrr);
    fprintf(fid, '\r\n\r\n');
    fprintf(fid, 'Calibration Images: ');
    fprintf(fid, imagefilestrr);
    fprintf(fid, '\r\n');
    fprintf(fid, 'Calibration Images Path: ');
    fprintf(fid, '\r\n');
    fprintf(fid, '%s', dosePath);
    fclose(fid);
   

end

