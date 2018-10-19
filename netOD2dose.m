function dosedir = netOD2dose(nODfiles, path, polyorder, p)
%convert netODFiles to Dose files

%create the dose file folder
idcs                    = strfind(path,'\');
parent                  = path(1:idcs(end-1));
dosedir                 = [parent, 'Dose\'];
[status, msg, msdID]    = mkdir(dosedir);

for n = 1:numel(nODfiles)
    netODfilename       = fullfile(path,nODfiles{n});
    netODArray          = imread(netODfilename);
    doseArray           = polyval(p, netODArray);
    dosefilename        = ['Dose_', nODfiles{n}];
    dosefullPath        = fullfile(dosedir, dosefilename);
     
    tiffwrite(doseArray, dosefullPath);
end

end


function tiffwrite(im, imname)
    
    if isa(im, 'uint16')
        tagstruct.BitsPerSample   = 16;
        tagstruct.SampleFormat    = Tiff.SampleFormat.UInt;
    elseif isa(im, 'single')
        tagstruct.BitsPerSample   = 32;
        tagstruct.SampleFormat    = Tiff.SampleFormat.IEEEFP;
    end
    
    t                         = Tiff(imname, 'w');
    tagstruct.ImageLength     = size(im, 1);
    tagstruct.ImageWidth      = size(im, 2);
    tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
    
    tagstruct.SamplesPerPixel = 1; 
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Compression     = Tiff.Compression.None; 
    setTag(t, tagstruct);
    t.write(im);
    t.close();
end