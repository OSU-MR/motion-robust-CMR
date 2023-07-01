function [ mag, phase_1, phase_2, phase_3] = saveDiccom4D_new( directory_name, info_dicom, data)
%READDICCOM4D Summary of this function goes here
%   Detailed explanation goes here

% Create a new folder for the DICOM
mkdir('DICOM')
cd('DICOM')
mkdir('MAG')
cd('MAG')
cd('..')

mkdir('X')
cd('X')
cd('..')

mkdir('Y')
cd('Y')
cd('..')

mkdir('Z')
cd('Z')
cd('..')

cur_dir = pwd;
cd(directory_name)

% Reshpae the data
data.mag = uint16((data.mag));
data.phase_1 = uint16((data.phase_1));
data.phase_2 = uint16((data.phase_2));
data.phase_3 = uint16((data.phase_3));

cd('MAG')
% info_dicom.WindowCenter= 26;
% info_dicom.WindowWidth= 5210;
info_dicom.SequenceName = 'fl3d1r2';  % Set for Magnitdue in 4D flow
info_dicom.SeriesDescription = '4DFLOW_nav_785_v2_GRAPPA_FB';
info_dicom.MediaStorageSOPInstanceUID = dicomuid;
info_dicom.SeriesNumber = 1;
info_dicom.StudyInstanceUID = dicomuid;
info_dicom.SeriesInstanceUID = dicomuid;
info_dicom.SOPInstanceUID = dicomuid;
info_dicom.ImageType = 'ORIGINAL\PRIMARY\M\ND';
saveImage(data.mag,info_dicom);
cd ..

cd 'X'
% info_dicom.WindowCenter= 2^11;
% info_dicom.WindowWidth= 2^12-1;
% info_dicom.LargestImagePixelValue= 2^12-1;
info_dicom.BitDepth = 12;
info_dicom.BitsStored = 12;
info_dicom.HighBit = 11;
info_dicom.LargestImagePixelValue=4094;
info_dicom.WindowCenter=26;
info_dicom.WindowWidth=4096;
info_dicom.RescaleIntercept= -4096;
info_dicom.RescaleSlope= 2;
info_dicom.RescaleType= 'US';
info_dicom.WindowCenterWidthExplanation= 'Algo1';
info_dicom.SequenceName = 'fl3d1_v350in';  % Set for Phase1
info_dicom.SeriesDescription = '4DFLOW_nav_785_v2_GRAPPA_FB_P';
info_dicom.SeriesNumber = 2;
info_dicom.SeriesInstanceUID = dicomuid;
info_dicom.SOPInstanceUID = dicomuid;
info_dicom.ImageType = 'DERIVED\PRIMARY\P\ND';
saveImage(data.phase_1,info_dicom);
cd ..

cd 'Y'
% info_dicom.WindowCenter= 2^11;
% info_dicom.WindowWidth= 2^12-1;
% info_dicom.BitDepth = 12;
% info_dicom.BitsStored = 12;
% info_dicom.LargestImagePixelValue= 2^12-1;
info_dicom.SequenceName = 'fl3d1_v350ap';  % Set for Phase2
info_dicom.SeriesDescription = '4DFLOW_nav_785_v2_GRAPPA_FB_P';
info_dicom.SeriesNumber = 3;
info_dicom.SeriesInstanceUID = dicomuid;
info_dicom.SOPInstanceUID = dicomuid;
info_dicom.ImageType = 'DERIVED\PRIMARY\P\ND';
saveImage(data.phase_2,info_dicom);
cd ..

cd 'Z'
% info_dicom.WindowCenter= 2^11;
% info_dicom.WindowWidth= 2^12-1;
% info_dicom.BitDepth = 12;
% info_dicom.BitsStored = 12;
% info_dicom.LargestImagePixelValue= 2^12-1;
info_dicom.SequenceName = 'fl3d1_v350fh';  % Set for Phase3
info_dicom.SeriesDescription = '4DFLOW_nav_785_v2_GRAPPA_FB_P';
info_dicom.SeriesNumber = 4;
info_dicom.SeriesInstanceUID = dicomuid;
info_dicom.SOPInstanceUID = dicomuid;
info_dicom.ImageType = 'DERIVED\PRIMARY\P\ND';
saveImage(data.phase_3,info_dicom);
cd ..

cd(cur_dir)


function saveImage(image,info_dicom)
    % Loop through all the indices
    fprintf(sprintf('\nSaving data   %0.1f%%%%',0.0))
    instanceID = 1;
    info_dicom.SliceLocation = 0;
    info_dicom.InstanceNumber = instanceID;
    
    for slice_ind = 1:size(image,3)
        info_dicom.ImagePositionPatient = [polyval([0.042442986126332   1.121072788336742],slice_ind),...
                                           polyval([-0.011314025097300  -1.311377868143300],slice_ind),...
                                           polyval([-0.040872799815748   2.036225159537864],slice_ind)];
        % Set the Trigger time to zero for each series of frames
        info_dicom.TriggerTime = 0;
        for frame_ind = 1:size(image,4)
            dicomwrite(image(:,:,slice_ind,frame_ind), ['MID000',num2str(instanceID)],...
                info_dicom);
            
            % Increment trigger time by one TR each frame
            info_dicom.TriggerTime = info_dicom.TriggerTime + info_dicom.RepetitionTime;
%             fprintf(sprintf('trigger tiem %f\n',info_dicom.TriggerTime));
            % Increment instance ID
            instanceID = instanceID +1;
            info_dicom.InstanceNumber = instanceID;
        end
        info_dicom.SliceLocation = info_dicom.SliceLocation +...
                                info_dicom.SliceThickness;
        info_dicom.ImagePositionPatient(3) = info_dicom.ImagePositionPatient(3) +1;
        
        % Display progress
        if length(sprintf('%0.1f%%%%',100*slice_ind/size(image,3)))>5
            fprintf(sprintf('\b\b\b\b\b'))
            fprintf(sprintf('%0.1f%%%%',100*slice_ind/size(image,3)))
        else
            fprintf(sprintf('\b\b\b\b'))
            fprintf(sprintf('%0.1f%%%%',100*slice_ind/size(image,3)))
        end
    end

end

end

