function saveDicom2Dt( x, dicom_info, fname )
%SAVEDICCOM2DT saves complex valued image sequence x, with dicom header
%dicom_info to filename fname
%   Detailed explanation goes here

image = uint16(2^12*x);

cur_dir = pwd;
mkdir(fname);
cd(fname);
tag = num2str(abs(round(100*randn(1))));
% Save Magnitude image
trigger = dicom_info.TriggerTime;
h = waitbar(0,'Saving Image Series');
for ind = 1:size(x,3)
    info = dicom_info;
    info.InstanceNumber = ind;
%     info.TriggerTime = trigger;
%     info.WindowCenter = 2047;
%     info.WindowWidth = 4095;
%     info.PixelSpacing = [2.34;3];
    dicomwrite(image(:,:,ind),['MID',tag,num2str(ind)],info,'CreateMode','Copy','WritePrivate',1)
%     dicomwrite(image(:,:,ind),['MID',tag,num2str(ind)],info);
% %     'RTFLOW_TEST.MR._.0024.0001.2017.10.19.16.33.20.295464.80174391.IMA'
    if ~mod(ind,10)
        waitbar(ind/size(x,3),h)
    end
    trigger = trigger+info.RepetitionTime;
end
close(h)
cd ..

cd(cur_dir)
end

