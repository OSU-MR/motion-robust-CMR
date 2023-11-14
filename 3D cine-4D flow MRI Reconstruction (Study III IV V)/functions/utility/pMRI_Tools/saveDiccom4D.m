function [ mag, phase_1, phase_2, phase_3] = saveDiccom4D( directory_name, dirs, data)
%READDICCOM4D Summary of this function goes here
%   Detailed explanation goes here

% Create a new folder for the DICOM
mkdir('DICOM')
cd('DICOM')
mkdir(dirs.mag)
cd(dirs.mag)
save_dir.mag = pwd;
cd('..')

mkdir(dirs.phase_1)
cd(dirs.phase_1)
save_dir.phase_1 = pwd;
cd('..')

mkdir(dirs.phase_2)
cd(dirs.phase_2)
save_dir.phase_2 = pwd;
cd('..')

mkdir(dirs.phase_3)
cd(dirs.phase_3)
save_dir.phase_3 = pwd;
cd('../..')

cur_dir = pwd;
cd(directory_name)

%% Copy the Dicomdir file
copyfile('DICOMDIR',[save_dir.mag,'/../DICOMDIR'])

% Reshpae the data
sizes = size(data.mag);
data.mag = uint16(reshape(data.mag,[sizes(1), sizes(2), sizes(3)*sizes(4)]));
data.phase_1 = uint16(reshape(data.phase_1,[sizes(1), sizes(2), sizes(3)*sizes(4)]));
data.phase_2 = uint16(reshape(data.phase_2,[sizes(1), sizes(2), sizes(3)*sizes(4)]));
data.phase_3 = uint16(reshape(data.phase_3,[sizes(1), sizes(2), sizes(3)*sizes(4)]));

cd([directory_name,'/',dirs.mag])
saveImage(data.mag,save_dir.mag);

cd([directory_name,'/',dirs.phase_1])
saveImage(data.phase_1,save_dir.phase_1);

cd([directory_name,'/',dirs.phase_2])
saveImage(data.phase_2,save_dir.phase_2);

cd([directory_name,'/',dirs.phase_3])
saveImage(data.phase_3,save_dir.phase_3);

cd(cur_dir)


function saveImage(image,save_dir)
    dir_current = pwd;
    % Loop through the Diccom images 
    direct = dir(pwd);

    % Get som info from the Dicom
    num_images = length(direct) - 2;

    %% Import the data into a 4D array 

    % loop through all images in the directory and import them in the Instance
    % number order
    fprintf(sprintf('\nSaving data   %0.1f%%%%',0.0))

    for ind = 3:length(direct)
        
        % Get dicominfo
        info = dicominfo(direct(ind).name);  
%         info = copy_info(info);
%         image_no_order(:,:,info.InstanceNumber) = double(dicomread(direct(ind).name));
        cd(save_dir)
        dicomwrite(image(:,:,info.InstanceNumber), direct(ind).name, info, 'CreateMode', 'copy');
%         dicomwrite(image(:,:,info.InstanceNumber), direct(ind).name, info);
        cd(dir_current)
        position(info.InstanceNumber,:) = info.ImagePositionPatient;
%         if ind ==20
%                 fprintf(sprintf('OR1 %f  ',info.ImageOrientationPatient(1:3)));
%                 fprintf('\n');
%                 fprintf(sprintf('OR2 %f  ',info.ImageOrientationPatient(4:6)));
%                 fprintf('\n');
%                 fprintf('\n');
                
                fprintf(sprintf('OR1 %f  ',info.ImagePositionPatient));
%                 fprintf('\n');
%                 fprintf(sprintf('OR2 %f  ',info.ImageOrientationPatient(4:6)));
                fprintf('\n');
                fprintf('\n');
%                 fprintf(sprintf('%s\n',info.ScanningSequence));
%                 fprintf(sprintf('%s\n',info.ScanOptions));
%                 fprintf(sprintf('%s\n',info.MRAcquisitionType));
%                 fprintf(sprintf('%s\n',info.SequenceName));
%                 fprintf(sprintf('%s\n',info.StationName));
%                 fprintf(sprintf('%s\n',info.StudyDescription));
%                 fprintf(sprintf('%s\n',info.SeriesDescription));
%                 fprintf(sprintf('%s\n',info.Manufacturer));
%                 info.SequenceVariant
%                 fprintf('\n');
                
%         end
%         % Display progress
%         if ~mod(ind,20)
%             if length(sprintf('%0.1f%%%%',100*ind/num_images))>5
%                 fprintf(sprintf('\b\b\b\b\b'))
%                 fprintf(sprintf('%0.1f%%%%',100*ind/num_images))
%             else
%                 fprintf(sprintf('\b\b\b\b'))
%                 fprintf(sprintf('%0.1f%%%%',100*ind/num_images))
%             end
%         end
    end
    fprintf(sprintf('\b\b\b\b\b'))
    fprintf(sprintf('%0.1f%%%%',100))
    

end

    function info_new = copy_info(info)
%         info_new.Filename = info.Filename;
%         info_new.FileModDate = info.FileModDate;
%         info_new.FileSize = info.FileSize;
        info_new.Format = info.Format;
        info_new.FormatVersion = info.FormatVersion;
        info_new.Width = info.Width;
        info_new.Height = info.Height;
        info_new.BitDepth = info.BitDepth;
        info_new.ColorType = info.ColorType;
%     info_new.FileMetaInformationGroupLength = info.FileMetaInformationGroupLength;
%     info_new.FileMetaInformationVersion = info.FileMetaInformationVersion;
%     info_new.MediaStorageSOPClassUID = info.MediaStorageSOPClassUID;
%     info_new.MediaStorageSOPInstanceUID = info.MediaStorageSOPInstanceUID;
%     info_new.TransferSyntaxUID = info.TransferSyntaxUID;
%     info_new.ImplementationClassUID = info.ImplementationClassUID;
%     info_new.ImplementationVersionName = info.ImplementationVersionName;
%     info_new.SpecificCharacterSet = info.SpecificCharacterSet;
        info_new.ImageType = info.ImageType;
%     info_new.InstanceCreationDate = info.InstanceCreationDate;
%     info_new.InstanceCreationTime = info.InstanceCreationTime;
        info_new.SOPClassUID = info.SOPClassUID;
        info_new.SOPInstanceUID = info.SOPInstanceUID;
%     info_new.StudyDate = info.StudyDate;
%     info_new.SeriesDate = info.SeriesDate;
%     info_new.AcquisitionDate = info.AcquisitionDate;
%     info_new.ContentDate = info.ContentDate;
%     info_new.StudyTime = info.StudyTime;
%     info_new.SeriesTime = info.SeriesTime;
%     info_new.AcquisitionTime = info.AcquisitionTime;
%     info_new.ContentTime = info.ContentTime;
        info_new.AccessionNumber = info.AccessionNumber;
        info_new.Modality = info.Modality;
    info_new.Manufacturer = info.Manufacturer;
%     info_new.InstitutionName = info.InstitutionName;
%     info_new.InstitutionAddress = info.InstitutionAddress;
%     info_new.ReferringPhysicianName = info.ReferringPhysicianName;
    info_new.StationName = info.StationName;
    info_new.StudyDescription = info.StudyDescription;
    info_new.SeriesDescription = info.SeriesDescription;
%     info_new.InstitutionalDepartmentName = info.InstitutionalDepartmentName;
%     info_new.OperatorName = info.OperatorName;
        info_new.ManufacturerModelName = info.ManufacturerModelName;
%     info_new.ReferencedImageSequence = info.ReferencedImageSequence;
%     info_new.PatientName = info.PatientName;
%     info_new.PatientID = info.PatientID;
%     info_new.PatientBirthDate = info.PatientBirthDate;
%     info_new.PatientSex = info.PatientSex;
%     info_new.PatientAge = info.PatientAge;
%     info_new.PatientSize = info.PatientSize;
%     info_new.PatientWeight = info.PatientWeight;
    info_new.BodyPartExamined = info.BodyPartExamined;
    info_new.ScanningSequence = info.ScanningSequence;
    info_new.SequenceVariant = info.SequenceVariant;
    info_new.ScanOptions = info.ScanOptions;
    info_new.MRAcquisitionType = info.MRAcquisitionType;
    info_new.SequenceName = info.SequenceName;
        info_new.AngioFlag = info.AngioFlag;
        info_new.SliceThickness = info.SliceThickness;
        info_new.RepetitionTime = info.RepetitionTime;
        info_new.EchoTime = info.EchoTime;
        info_new.NumberOfAverages = info.NumberOfAverages;
        info_new.ImagingFrequency = info.ImagingFrequency;
%     info_new.ImagedNucleus = info.ImagedNucleus;
        info_new.EchoNumber = info.EchoNumber;
        info_new.MagneticFieldStrength = info.MagneticFieldStrength;
        info_new.NumberOfPhaseEncodingSteps = info.NumberOfPhaseEncodingSteps;
        info_new.EchoTrainLength = info.EchoTrainLength;
%     info_new.PercentSampling = info.PercentSampling;
%     info_new.PercentPhaseFieldOfView = info.PercentPhaseFieldOfView;
        info_new.PixelBandwidth = info.PixelBandwidth;
%     info_new.DeviceSerialNumber = info.DeviceSerialNumber;
%     info_new.SoftwareVersion = info.SoftwareVersion;
%     info_new.ProtocolName = info.ProtocolName;
        info_new.TriggerTime = info.TriggerTime;
        info_new.NominalInterval = info.NominalInterval;
        info_new.CardiacNumberOfImages = info.CardiacNumberOfImages;
        info_new.TransmitCoilName = info.TransmitCoilName;
        info_new.AcquisitionMatrix = info.AcquisitionMatrix;
        info_new.InPlanePhaseEncodingDirection = info.InPlanePhaseEncodingDirection;
        info_new.FlipAngle = info.FlipAngle;
        info_new.VariableFlipAngleFlag = info.VariableFlipAngleFlag;
%     info_new.SAR = info.SAR;
%     info_new.dBdt = info.dBdt;
%     info_new.PatientPosition = info.PatientPosition;
%         info_new.Private_0019_10xx_Creator = info.Private_0019_10xx_Creator;
%         info_new.Private_0019_1008 = info.Private_0019_1008;
%         info_new.Private_0019_1009 = info.Private_0019_1009;
%         info_new.Private_0019_100b = info.Private_0019_100b;
%         info_new.Private_0019_100f = info.Private_0019_100f;
%         info_new.Private_0019_1011 = info.Private_0019_1011;
%         info_new.Private_0019_1012 = info.Private_0019_1012;
%         info_new.Private_0019_1013 = info.Private_0019_1013;
%         info_new.Private_0019_1014 = info.Private_0019_1014;
%         info_new.Private_0019_1015 = info.Private_0019_1015;
%         info_new.Private_0019_1017 = info.Private_0019_1017;
%         info_new.Private_0019_1018 = info.Private_0019_1018;
        info_new.StudyInstanceUID = info.StudyInstanceUID;
        info_new.SeriesInstanceUID = info.SeriesInstanceUID;
        info_new.StudyID = info.StudyID;
        info_new.SeriesNumber = info.SeriesNumber;
        info_new.AcquisitionNumber = info.AcquisitionNumber;
        info_new.InstanceNumber = info.InstanceNumber;
        info_new.ImagePositionPatient = info.ImagePositionPatient;
        info_new.ImageOrientationPatient = info.ImageOrientationPatient;
        
        info_new.FrameOfReferenceUID = info.FrameOfReferenceUID;
        info_new.PositionReferenceIndicator = info.PositionReferenceIndicator;
        info_new.SliceLocation = info.SliceLocation;
        info_new.SamplesPerPixel = info.SamplesPerPixel;
        info_new.PhotometricInterpretation = info.PhotometricInterpretation;
        info_new.Rows = info.Rows;
        info_new.Columns = info.Columns;
        info_new.PixelSpacing = info.PixelSpacing;
        info_new.BitsAllocated = info.BitsAllocated;
        info_new.BitsStored = info.BitsStored;
        info_new.HighBit = info.HighBit;
        info_new.PixelRepresentation = info.PixelRepresentation;
        info_new.SmallestImagePixelValue = info.SmallestImagePixelValue;
        info_new.LargestImagePixelValue = info.LargestImagePixelValue;
        info_new.WindowCenter = info.WindowCenter;
        info_new.WindowWidth = info.WindowWidth;
        info_new.WindowCenterWidthExplanation = info.WindowCenterWidthExplanation;
%         info_new.Private_0029_10xx_Creator = info.Private_0029_10xx_Creator;
%         info_new.Private_0029_11xx_Creator = info.Private_0029_11xx_Creator;
%         info_new.Private_0029_1008 = info.Private_0029_1008;
%         info_new.Private_0029_1009 = info.Private_0029_1009;
%         info_new.Private_0029_1010 = info.Private_0029_1010;
%         info_new.Private_0029_1018 = info.Private_0029_1018;
%         info_new.Private_0029_1019 = info.Private_0029_1019;
%         info_new.Private_0029_1020 = info.Private_0029_1020;
%         info_new.Private_0029_1160 = info.Private_0029_1160;
%         info_new.RequestedProcedureDescription = info.RequestedProcedureDescription;
%         info_new.PerformedProcedureStepStartDate = info.PerformedProcedureStepStartDate;
%         info_new.PerformedProcedureStepStartTime = info.PerformedProcedureStepStartTime;
%         info_new.PerformedProcedureStepID = info.PerformedProcedureStepID;
%         info_new.PerformedProcedureStepDescription = info.PerformedProcedureStepDescription;
%         info_new.Private_0051_10xx_Creator = info.Private_0051_10xx_Creator;
%         info_new.Private_0051_1008 = info.Private_0051_1008;
%         info_new.Private_0051_1009 = info.Private_0051_1009;
%         info_new.Private_0051_100a = info.Private_0051_100a;
%         info_new.Private_0051_100b = info.Private_0051_100b;
%         info_new.Private_0051_100c = info.Private_0051_100c;
%         info_new.Private_0051_100d = info.Private_0051_100d;
%         info_new.Private_0051_100e = info.Private_0051_100e;
%         info_new.Private_0051_100f = info.Private_0051_100f;
%         info_new.Private_0051_1011 = info.Private_0051_1011;
%         info_new.Private_0051_1012 = info.Private_0051_1012;
%         info_new.Private_0051_1013 = info.Private_0051_1013;
%         info_new.Private_0051_1016 = info.Private_0051_1016;
%         info_new.Private_0051_1017 = info.Private_0051_1017;
%         info_new.Private_0051_1018 = info.Private_0051_1018;
%         info_new.Private_0051_1019 = info.Private_0051_1019;
%         fprintf(sprintf('%s\n',info.SequenceVariant));
%         info.SequenceVariant
%         fprintf(sprintf('%s\n',info_new.ImagePositionPatient));
%         fprintf(sprintf('%s\n',info_new.ImageOrientationPatient));
%         fprintf(sprintf('%s\n',info_new.SliceLocation));
%         display('');
    end
end

