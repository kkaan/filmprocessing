function [Cone_Data, MLC_Data, LinacName, Meas_Date, Norm_Data, Data_Dir] = filmexcelextract(filename)
%EXELEXTRACT Extracts required variables from excel data sheet
%   Parse the excel sheet with film dose descriptions
       
   
        [Meas_Date, ~, ~]   = xlsread(filename, "sheet1", "Meas_Date");
        [~, Data_Dir, ~]    = xlsread(filename, "sheet1", "Data_Dir");
        [~, LinacName, ~]   = xlsread(filename, "sheet1", "LinacName");
        [~, ~, MLC_Data]    = xlsread(filename, "sheet1", "MLC_Data");
        [~, ~, Norm_Data]   = xlsread(filename, "sheet1", "Norm_Data");
        [~, ~, Cone_Data]   = xlsread(filename, "sheet1", "Cone_Data");
      
    
end

