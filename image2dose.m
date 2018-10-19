

[doseFiles, dosePath, ~] = uigetfile({'*.tiff', 'Tiff Files'; '*.*', 'All Files'}, ...
                                                'Select Calibration Files', 'MultiSelect', 'On');
                                      

                                            
                                        
coeffstrr = sprintf('%.5f ' , app.p);
dosestrr = sprintf('%.2f ', doses);
imagefilestrr = sprintf('%s, ', doseFiles{:});

filename = fullfile(dosePath, 'caliFileName.txt');
fid = fopen(filename, 'w');
fprintf(fid, 'Polynomial Coefficients: ');
fprintf(fid, coeffstrr);
fprintf(fid, '\r\n');
fprintf(fid, 'Doses: ');
fprintf(fid, dosestrr);
fprintf(fid, 'Calibration Images: ');
fprintf(fid, imagefilestrr);
fprintf(fid, '\r\n');
fprintf(fid, '%s', dosePath);
fclose(fid);

                                            