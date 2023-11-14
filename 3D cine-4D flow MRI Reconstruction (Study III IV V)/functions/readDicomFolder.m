% Read all files from the dicom folder
% clear all,
function [a_0, t_info, N_S, fnames] = readDicomFolder(pname)

%[fname,pname] = uigetfile('*.*','Select one file among all filtes that you want to read');
fnames = dir(pname);
l_name = length(fnames);
j = 0;
for i=1:l_name % find all the files
    if (fnames(i).isdir==0)&(~strcmp(fnames(i).name,'DICOMDIR'))&(~strcmp(fnames(i).name,'dicomdir'))&(~strcmp(fnames(i).name,'Dicomdir')) 
        n_j(j+1) = i;
        j = j+1;
    end
end

l_names = length(n_j); % file names


% h = waitbar(0,'Loading Dicom files. Please wait...', 'name', 'Loading');
t0 = dicomread( sprintf('%s%s',pname,fnames(n_j(2)).name) );
s = size(t0);
a_0 = zeros(s(1),s(2),l_names); % initialize the original data matrix


k = 1 ;
s_order = 0;
for i=1:l_names
    t_info(i) = dicominfo( sprintf( '%s%s',pname,fnames(n_j(i)).name ) ); 
    temp = dicomread( sprintf( '%s%s',pname,fnames(n_j(i)).name ) );
    %a_mean(i) = mean(temp( : )); % Find the mean of each image
    s_t =  size(temp);
    if (s(1) == s_t(1)) 
        a_0(:,:,i) = temp  ; % assign the pixel value to a 3D matrix
    else
        s_order(k) = i; % Deal with transposed matrix
        a_0(:,:,i) = temp' ;
        k = k + 1;
    end
    if ischar(t_info(i).InstanceNumber)
        Instance_Number(i) = str2num(t_info(i).InstanceNumber);
    else
        Instance_Number(i) = t_info(i).InstanceNumber;
    end
%     waitbar(i/l_names);
end

% close(h)
s = size(a_0);
[s_In, s_Order ] = sort(Instance_Number); % Reorder them accordig to instance number
a_0 = a_0(:,:,s_Order); % Sort a_0
t_info = t_info(s_Order); % sort t_info
%size(s_order),
%size(s_Order),
if s_order==0, % without transpose slices
else % with transpose slices
    for j=1:length(s_order),s_order(j) = find(s_Order == s_order(j));end % find the new transposed images
end
for i=1:s(3),sl(i) = t_info(i).SliceLocation;end
sl_o = sort(sl); % ordered
sl_d = diff(sl_o);% difference
N_S = sum(abs(sl_d)>0.01)+1; % Number of slices
sl(s_Order);



%s = size(a_0);
% clear temp; clear h, clear t0; clear s, clear l_name; clear i, clear j; clear n_j; 
% [filename, pathname] = uiputfile('*.mat', 'Save Workspace as');
% save(sprintf('%s%s',pathname, filename))
