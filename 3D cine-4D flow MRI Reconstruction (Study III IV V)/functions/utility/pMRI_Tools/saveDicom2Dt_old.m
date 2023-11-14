function [ output_args ] = saveDicom2Dt( x, dicom_info, fname )
%SAVEDICCOM2DT saves complex valued image sequence x, with dicom header
%dicom_info to filename fname
%   Detailed explanation goes here

mag = uint16(2^12*abs(x)/max(abs(x(:))));
phase = uint16(2^12*norm_velocity(angle(x)));

cur_dir = pwd;
mkdir(fname);
cd(fname);

% Save Magnitude image
mkdir([fname,'_mag']);
cd([fname,'_mag']);
trigger = dicom_info.TriggerTime;
h = waitbar(0,'Saving magnitude image');
for ind = 1:size(x,3)
    info = dicom_info;
    info.InstanceNumber = ind;
    info.SeriesNumber = 1;
    info.TriggerTime = trigger;
    trigger
    dicomwrite(mag(:,:,ind),['MID00',num2str(ind)],info);
    if ~mod(ind,10)
        waitbar(ind/size(x,3),h)
    end
    trigger = trigger+info.RepetitionTime;
end
close(h)
cd ..

% Save Phase image
mkdir([fname,'_phase']);
cd([fname,'_phase']);
h = waitbar(0,'Saving phase image');
trigger = dicom_info.TriggerTime;
for ind = 1:size(x,3)
    info = dicom_info;
    info.InstanceNumber = ind;
    info.SeriesNumber = 2;
    info.TriggerTime = trigger;
    trigger
    dicomwrite(phase(:,:,ind),['MID00',num2str(ind+size(x,3))],info);
    if ~mod(ind,10)
        waitbar(ind/size(x,3),h)
    end
    trigger = trigger+info.RepetitionTime;
end
close(h)
cd(cur_dir)
end

