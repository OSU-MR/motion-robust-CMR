function [ data, dirs ] = readDiccom4D( directory_name )
%READDICCOM4D Summary of this function goes here
%   Detailed explanation goes here

cur_dir = pwd;
directory = dir(directory_name);
cd(directory_name)

data.mag = [];
data.phase_1 = [];
data.phase_2 = [];
data.phase_3 = [];
phase_count = 0;

dirs.mag = [];
dirs.phase_1 = [];
dirs.phase_2 = [];
dirs.phase_3 = [];

% Read all directories
% warning('importing every other directory, change if not desired');
for dir_ind = 3:length(directory)
    folder = dir(directory(dir_ind).name);
    
    % Make sure we are looking a directory
    if directory(dir_ind).isdir
        cd(directory(dir_ind).name);
        dicom_info = dicominfo(folder(3).name); % Dicom info for the first frame

        % if magnitude image
        if isempty(strfind(dicom_info.Private_0051_1016,'P/'))
            fprintf('\nmagnitude image detected')
            data.mag = readImage();
            dirs.mag = directory(dir_ind).name;
        % If Phase Image
        else
            fprintf('\nphase image detected')
            if phase_count == 0
                data.phase_1 = readImage();
                dirs.phase_1 = directory(dir_ind).name;
            elseif phase_count == 1
                data.phase_2 = readImage();
                dirs.phase_2 = directory(dir_ind).name;
            elseif phase_count == 2                
                data.phase_3 = readImage();
                dirs.phase_3 = directory(dir_ind).name;
            end
            phase_count = phase_count +1;

        end
        cd ..
    end
end

%               switch lower(dicom_info.Private_0051_1014(end-1:end))
%               case 'rl'
%                 fprintf('\nrl phase image detected')
%                     directory(dir_ind).name
%                     dicom_info.Private_0051_1016
% %                 phase_rl = readImage(directory(dir_ind).name);
%               case 'ap'
%                 fprintf('\nap image detected')
%                     directory(dir_ind).name
%                     dicom_info.Private_0051_1016
% %                 phase_ap = readImage(directory(dir_ind).name);
%               case 'gh'
%                 fprintf('\nthrough plane phase image detected')
%                     directory(dir_ind).name
%                     dicom_info.Private_0051_1016
% %                 phase_through = readImage(directory(dir_ind).name);
%               otherwise
%                 disp('Unknown Phase Direction')
%               end    
% num_folders = length(directory)-3;

% image1 = readImage([directory_name,'/SER_',sprintf('%d',num_folders-30)]); %44
% image2 = readImage([directory_name,'/SER_',sprintf('%d',num_folders-31)]); %43
% image3 = readImage([directory_name,'/SER_',sprintf('%d',num_folders-32)]); %42
% image4 = readImage([directory_name,'/SER_',sprintf('%d',num_folders-33)]); %41
% for ii = 7:-1:0
%     image(:,:,:,:,ii) = readImage([directory_name,'/SER_',sprintf('%d',num_folders-ii)]);
% end

cd(cur_dir)
function image = readImage()
    % Loop through the Diccom images 
    direct = dir(pwd);

    % Get som info from the Dicom
    info = dicominfo(direct(3).name); % Dicom info for the first frame
    num_images = length(direct) - 2;
    num_frames = info.CardiacNumberOfImages;
    num_slices = num_images/num_frames;

    image_no_order = zeros(info.Rows, info.Columns, num_images);
    image = zeros(info.Rows, info.Columns, num_slices ,num_frames);
    %% Import the data into a 4D array 

    % loop through all images in the directory and import them in the Instance
    % number order
    fprintf(sprintf('\nImporting data   %0.1f%%%%',0.0))
    for ind = 3:length(direct);
        
        % Get dicominfo
        info = dicominfo(direct(ind).name);  

        % Read dicom file and store it in a 3D array
        image_no_order(:,:,info.InstanceNumber) = double(dicomread(direct(ind).name));
        if ~mod(ind,20)
            if length(sprintf('%0.1f%%%%',100*ind/num_images))>5
                fprintf(sprintf('\b\b\b\b\b'))
                fprintf(sprintf('%0.1f%%%%',100*ind/num_images))
            else
                fprintf(sprintf('\b\b\b\b'))
                fprintf(sprintf('%0.1f%%%%',100*ind/num_images))
            end
        end
    end

    % Loop through the 3D no order array and place it in a 4D array
    ind = 1;
    fprintf(sprintf('\nReordering data...'))
    for frame_ind = 1:num_frames
        for slice_ind = 1:num_slices 
            image(:,:,slice_ind, frame_ind) = image_no_order(:,:,ind);
            ind = ind+1;
        end
    end
    fprintf(sprintf('done!\n'))
    
end

end


