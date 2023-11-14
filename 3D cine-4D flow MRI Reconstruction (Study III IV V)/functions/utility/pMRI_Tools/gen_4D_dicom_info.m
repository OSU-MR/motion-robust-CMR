function info = gen_4D_dicom_info(image,parameters)

%% Things to Set Per image
info.Width= size(image,2);
info.Columns= size(image,2);
info.Height= size(image,1);
info.Rows= size(image,1);

info.CardiacNumberOfImages= size(image,4);
info.RepetitionTime= parameters.TR;
info.EchoTime= parameters.TE;
info.PixelSpacing= parameters.pixelSpacing;
info.SliceThickness= parameters.sliceThickness;

info.SequenceName = 'fl3d1r2';  % Set for Magnitdue in 4D flow
% info.SequenceName = 'fl3d1_v350in';  % Set for Phase1
% info.SequenceName = 'fl3d1_v350ap';  % Set for Phase2
% info.SequenceName = 'fl3d1_v350fh';  % Set for Phase3
%% no sure how to seet these in general
info.ImagingFrequency= 63.622660000000003;
info.AcquisitionMatrix= [0,128,128,0];
% info.ImagePositionPatient= 1e2*[1.418173691221100,-1.390576043824400,1.750115560827600];
info.ImagePositionPatient = [0,0,0];
info.ImageOrientationPatient= [0.089461127746120, 0.979892418630680, -0.178346725600700, -0.701146099952100, -0.065217347244900, -0.710028762896500];
% info.ImageOrientationPatient= [1, 0, 0, 0, -1, 0];


%% Things to set per image frame and per slice
info.InstanceNumber= 1;   % start at 1 and increment for each image
info.SliceLocation= 0;      % Slice location, increment by slice thickness when saving new slices     
info.TriggerTime= 0;      % increment by repitition time from first image

%% Not Really important to set
info.FlipAngle= 10;
info.NumberOfAverages= 1;
info.PixelBandwidth= 505;
info.NominalInterval= 1025;
info.AcquisitionNumber= 1;
info.SamplesPerPixel= 1;
info.NumberOfPhaseEncodingSteps= 120;
info.SeriesNumber= 1;

%% Window Parameters may be important for phase
info.BitsAllocated= 16;
info.BitsStored= 16;
info.HighBit= 15;
info.PixelRepresentation= 0;
info.SmallestImagePixelValue= 0;
info.LargestImagePixelValue= 2^16-1;
info.WindowCenter= 162;
info.WindowWidth= 414;
info.WindowCenterWidthExplanation= 'Algo1';

info.BodyPartExamined = 'HEART';
info.ScanningSequence = 'GR';
info.ScanOptions = 'CT';
info.MRAcquisitionType = '3D';
info.StationName = 'MRC25383';
info.StudyDescription = 'OSU CMR LAB^Rizwan';
info.SeriesDescription = '4DFLOW_nav_785_v2_GRAPPA_FB_P';
info.Manufacturer = 'SIEMENS';
info.SequenceVariant = 'SK\SP\OSP';

%% Not really sure what these do
info.Format= 'DICOM';
info.FormatVersion= 3;
% info.BitDepth= 16;
info.ColorType= 'grayscale';
info.ImageType= 'ORIGINAL\PRIMARY\M\ND';
info.SOPClassUID= '1.2.840.10008.5.1.4.1.1.4';
info.SOPInstanceUID= '1.3.6.1.4.1.9590.100.1.2.7799186612403712300739662522884567978';
% info.AccessionNumber= '2577648E';
info.Modality= 'MR';
info.ManufacturerModelName= 'Avanto';
info.AngioFlag= 'N';
info.EchoNumber= 1;
info.MagneticFieldStrength= 1.500000000000000;
info.EchoTrainLength= 1;
info.TransmitCoilName= 'Body';
info.InPlanePhaseEncodingDirection= 'ROW';
info.VariableFlipAngleFlag= 'N';
info.StudyInstanceUID= '1.2.840.114350.2.172.2.798268.2.293062833.1';
info.SeriesInstanceUID= '1.3.12.2.1107.5.2.30.59096.2017071110365795484809777.0.0.0';
info.StudyID= '10';
info.FrameOfReferenceUID= '1.3.12.2.1107.5.2.30.59096.1.20170711095331546.0.0.0';
info.PositionReferenceIndicator= '';
info.PhotometricInterpretation= 'MONOCHROME2';



% %% Things to Set Per image
% info.Width= 224;
% info.Height= 320;
% info.Rows= 320;
% info.Columns= 224;
% info.CardiacNumberOfImages= 24;
% info.RepetitionTime= 35.200000000000003;
% info.EchoTime= 2.180000000000000;
% info.ImagingFrequency= 63.622660000000003;
% info.AcquisitionMatrix= [0,128,128,0];
% info.ImagePositionPatient= 1e2*[1.418173691221100,-1.390576043824400,1.750115560827600];
% info.ImageOrientationPatient= [0.089461127746120, 0.979892418630680, -0.178346725600700, -0.701146099952100, -0.065217347244900, -0.710028762896500];
% info.PixelSpacing= [0,0];
% info.SliceThickness= 6;
% 
% %% Things to set per image frame and per slice
% info.InstanceNumber= 1;   % start at 1 and increment for each image
% info.SliceLocation= 0;      % Slice location, increment by slice thickness when saving new slices     
% info.TriggerTime= 0;      % increment by repitition time from first image
% 
% %% Not Really important to set
% info.FlipAngle= 10;
% info.NumberOfAverages= 1;
% info.PixelBandwidth= 505;
% info.NominalInterval= 1025;
% info.AcquisitionNumber= 1;
% info.SamplesPerPixel= 1;
% info.NumberOfPhaseEncodingSteps= 120;
% info.SeriesNumber= 1;
% 
% %% Window Parameters may be important for phase
% info.BitsAllocated= 16;
% info.BitsStored= 16;
% info.HighBit= 15;
% info.PixelRepresentation= 0;
% info.SmallestImagePixelValue= 0;
% info.LargestImagePixelValue= 0;
% info.WindowCenter= 162;
% info.WindowWidth= 414;
% info.WindowCenterWidthExplanation= 'Algo1';
% 
% %% Not really sure what these do
% info.Format= 'DICOM';
% info.FormatVersion= 3;
% info.BitDepth= 16;
% info.ColorType= 'grayscale';
% info.ImageType= 'ORIGINAL\PRIMARY\M\ND';
% info.SOPClassUID= '1.2.840.10008.5.1.4.1.1.4';
% info.SOPInstanceUID= '1.3.6.1.4.1.9590.100.1.2.7799186612403712300739662522884567978';
% info.AccessionNumber= '2577648E';
% info.Modality= 'MR';
% info.ManufacturerModelName= 'Avanto';
% info.AngioFlag= 'N';
% info.EchoNumber= 1;
% info.MagneticFieldStrength= 1.500000000000000;
% info.EchoTrainLength= 1;
% info.TransmitCoilName= 'Body';
% info.InPlanePhaseEncodingDirection= 'ROW';
% info.VariableFlipAngleFlag= 'N';
% info.StudyInstanceUID= '1.2.840.114350.2.172.2.798268.2.293062833.1';
% info.SeriesInstanceUID= '1.3.12.2.1107.5.2.30.59096.2017071110365795484809777.0.0.0';
% info.StudyID= '293062833';
% info.FrameOfReferenceUID= '1.3.12.2.1107.5.2.30.59096.1.20170711095331546.0.0.0';
% info.PositionReferenceIndicator= '';
% info.PhotometricInterpretation= 'MONOCHROME2';
% 
