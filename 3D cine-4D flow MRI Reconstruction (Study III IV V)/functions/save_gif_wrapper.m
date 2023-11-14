function [] = save_gif_wrapper(obj, saveDir, saveName, fps, clip)
% Individual view cines
% -------------------------------------------------
obj_SAG = obj; obj_SAG.permuteImage([1,2,3,4]);
if ~exist([saveDir,'\sag\'],'dir')
    mkdir([saveDir,'\sag\']);
end
if ~exist([saveDir,'\cor\'],'dir')
    mkdir([saveDir,'\cor\']);
end
if ~exist([saveDir,'\tra\'],'dir')
    mkdir([saveDir,'\tra\']);
end

% SAG
% ****************************************************************
pattern = [1,2,3,4];
cB = permute(obj.outputs.xHatb, pattern);
cX = permute(obj.outputs.xHatx, pattern);
cY = permute(obj.outputs.xHaty, pattern);
cZ = permute(obj.outputs.xHatz, pattern);
tX = permute(obj.outputs.thetaX, pattern);
tY = permute(obj.outputs.thetaY, pattern);
tZ = permute(obj.outputs.thetaZ, pattern);
vX = permute(obj.outputs.vX, pattern);
vY = permute(obj.outputs.vY, pattern);
vZ = permute(obj.outputs.vZ, pattern);
% Normalize velocity
tX = norm_velocity(tX);
tY = norm_velocity(tY);
tZ = norm_velocity(tZ);
% Create gifs
create_GIF(abs(cB),[saveDir,'\sag\',saveName,'_sag_yB'], fps, clip);
create_GIF(abs(cX),[saveDir,'\sag\',saveName,'_sag_yX'], fps, clip);
create_GIF(abs(cY),[saveDir,'\sag\',saveName,'_sag_yY'], fps, clip);
create_GIF(abs(cZ),[saveDir,'\sag\',saveName,'_sag_yZ'], fps, clip);
create_GIF(tX,[saveDir,'\sag\',saveName,'_sag_thetaX'], fps, 1);
create_GIF(tY,[saveDir,'\sag\',saveName,'_sag_thetaY'], fps, 1);
create_GIF(tZ,[saveDir,'\sag\',saveName,'_sag_thetaZ'], fps, 1);
create_GIF(vX,[saveDir,'\sag\',saveName,'_sag_vX'], fps, 1);
create_GIF(vY,[saveDir,'\sag\',saveName,'_sag_vY'], fps, 1);
create_GIF(vZ,[saveDir,'\sag\',saveName,'_sag_vZ'], fps, 1);
% COR
% ****************************************************************
pattern = [1,3,2,4];
cB = permute(obj.outputs.xHatb, pattern);
cX = permute(obj.outputs.xHatx, pattern);
cY = permute(obj.outputs.xHaty, pattern);
cZ = permute(obj.outputs.xHatz, pattern);
tX = permute(obj.outputs.thetaX, pattern);
tY = permute(obj.outputs.thetaY, pattern);
tZ = permute(obj.outputs.thetaZ, pattern);
vX = permute(obj.outputs.vX, pattern);
vY = permute(obj.outputs.vY, pattern);
vZ = permute(obj.outputs.vZ, pattern);
% Normalize velocity
tX = norm_velocity(tX);
tY = norm_velocity(tY);
tZ = norm_velocity(tZ);
% Create gifs
create_GIF(abs(cB),[saveDir,'\cor\',saveName,'_cor_yB'], fps, clip);
create_GIF(abs(cX),[saveDir,'\cor\',saveName,'_cor_yX'], fps, clip);
create_GIF(abs(cY),[saveDir,'\cor\',saveName,'_cor_yY'], fps, clip);
create_GIF(abs(cZ),[saveDir,'\cor\',saveName,'_cor_yZ'], fps, clip);
create_GIF(tX,[saveDir,'\cor\',saveName,'_cor_thetaX'], fps, 1);
create_GIF(tY,[saveDir,'\cor\',saveName,'_cor_thetaY'], fps, 1);
create_GIF(tZ,[saveDir,'\cor\',saveName,'_cor_thetaZ'], fps, 1);
create_GIF(vX,[saveDir,'\cor\',saveName,'_cor_vX'], fps, 1);
create_GIF(vY,[saveDir,'\cor\',saveName,'_cor_vY'], fps, 1);
create_GIF(vZ,[saveDir,'\cor\',saveName,'_cor_vZ'], fps, 1);
% TRA
% ****************************************************************
pattern = [2,3,1,4];
cB = permute(obj.outputs.xHatb, pattern);
cX = permute(obj.outputs.xHatx, pattern);
cY = permute(obj.outputs.xHaty, pattern);
cZ = permute(obj.outputs.xHatz, pattern);
tX = permute(obj.outputs.thetaX, pattern);
tY = permute(obj.outputs.thetaY, pattern);
tZ = permute(obj.outputs.thetaZ, pattern);
vX = permute(obj.outputs.vX, pattern);
vY = permute(obj.outputs.vY, pattern);
vZ = permute(obj.outputs.vZ, pattern);
% Normalize velocity
tX = norm_velocity(tX);
tY = norm_velocity(tY);
tZ = norm_velocity(tZ);
% Create gifs
create_GIF(abs(cB),[saveDir,'\tra\',saveName,'_tra_yB'], fps, clip);
create_GIF(abs(cX),[saveDir,'\tra\',saveName,'_tra_yX'], fps, clip);
create_GIF(abs(cY),[saveDir,'\tra\',saveName,'_tra_yY'], fps, clip);
create_GIF(abs(cZ),[saveDir,'\tra\',saveName,'_tra_yZ'], fps, clip);
create_GIF(tX,[saveDir,'\tra\',saveName,'_tra_thetaX'], fps, 1);
create_GIF(tY,[saveDir,'\tra\',saveName,'_tra_thetaY'], fps, 1);
create_GIF(tZ,[saveDir,'\tra\',saveName,'_tra_thetaZ'], fps, 1);
create_GIF(vX,[saveDir,'\tra\',saveName,'_tra_vX'], fps, 1);
create_GIF(vY,[saveDir,'\tra\',saveName,'_tra_vY'], fps, 1);
create_GIF(vZ,[saveDir,'\tra\',saveName,'_tra_vZ'], fps, 1);
% ---------------------------------------------------
% Logical axis panel cines
% ---------------------------------------------------
% Sagittal 
sagB = permute(obj.outputs.xHatb, [1,2,3,4]);
sagX = permute(obj.outputs.thetaX, [1,2,3,4]);
sagY = permute(obj.outputs.thetaY, [1,2,3,4]);
sagZ = permute(obj.outputs.thetaZ, [1,2,3,4]);
% Coronal
corB = permute(obj.outputs.xHatb, [1,3,2,4]);
corX = permute(obj.outputs.thetaX, [1,3,2,4]);
corY = permute(obj.outputs.thetaY, [1,3,2,4]);
corZ = permute(obj.outputs.thetaZ, [1,3,2,4]);
% Transverse
traB = permute(obj.outputs.xHatb, [2,3,1,4]);
traX = permute(obj.outputs.thetaX, [2,3,1,4]);
traY = permute(obj.outputs.thetaY, [2,3,1,4]);
traZ = permute(obj.outputs.thetaZ, [2,3,1,4]);
% Choose central slice
sagB = squeeze(sagB(:,:,round(end/2),:));
sagX = squeeze(sagX(:,:,round(end/2),:));
sagY = squeeze(sagY(:,:,round(end/2),:));
sagZ = squeeze(sagZ(:,:,round(end/2),:));
corB = squeeze(corB(:,:,round(end/2),:));
corX = squeeze(corX(:,:,round(end/2),:));
corY = squeeze(corY(:,:,round(end/2),:));
corZ = squeeze(corZ(:,:,round(end/2),:));
traB = squeeze(traB(:,:,round(end/2),:));
traX = squeeze(traX(:,:,round(end/2),:));
traY = squeeze(traY(:,:,round(end/2),:));
traZ = squeeze(traZ(:,:,round(end/2),:));
% Take magnitude
sagB = abs(sagB);
corB = abs(corB);
traB = abs(traB);
% Clip magnitude
sagB(sagB > (1/clip)*max(sagB(:))) = (1/clip)*max(sagB(:));
corB(corB > (1/clip)*max(corB(:))) = (1/clip)*max(corB(:));
traB(traB > (1/clip)*max(traB(:))) = (1/clip)*max(traB(:));
% Rescale
sagB = rescale(sagB, 0, 1);
sagX = rescale(sagX, 0, 1);
sagY = rescale(sagY, 0, 1);
sagZ = rescale(sagZ, 0, 1);
corB = rescale(corB, 0, 1);
corX = rescale(corX, 0, 1);
corY = rescale(corY, 0, 1);
corZ = rescale(corZ, 0, 1);
traB = rescale(traB, 0, 1);
traX = rescale(traX, 0, 1);
traY = rescale(traY, 0, 1);
traZ = rescale(traZ, 0, 1);
% Define panels
panel = zeros([size(sagB,1)+size(corB,1)+size(traB,1), 4*max([size(sagB,2),size(corB,2),size(traB,2)]), size(sagB,3)]);
dimRows = [size(sagB,1), size(corB,1), size(traB,1)];
dimCols = [size(sagB,2), size(corB,2), size(traB,2)];
mx = max(dimCols);
center = [floor(mx-(mx/2))+1, floor(2*mx-(mx/2))+1, floor(3*mx-(mx/2))+1, floor(4*mx-(mx/2))+1];
% Build panels
panel(1:dimRows(1), center(1)-ceil(dimCols(1)/2)+1:center(1)+floor(dimCols(1)/2), :) = sagB;
panel(1:dimRows(1), center(2)-ceil(dimCols(1)/2)+1:center(2)+floor(dimCols(1)/2), :) = sagX;
panel(1:dimRows(1), center(3)-ceil(dimCols(1)/2)+1:center(3)+floor(dimCols(1)/2), :) = sagY;
panel(1:dimRows(1), center(4)-ceil(dimCols(1)/2)+1:center(4)+floor(dimCols(1)/2), :) = sagZ;
%
panel(dimRows(1)+1:dimRows(1)+dimRows(2), center(1)-ceil(dimCols(2)/2)+1:center(1)+floor(dimCols(2)/2), :) = corB;
panel(dimRows(1)+1:dimRows(1)+dimRows(2), center(2)-ceil(dimCols(2)/2)+1:center(2)+floor(dimCols(2)/2), :) = corX;
panel(dimRows(1)+1:dimRows(1)+dimRows(2), center(3)-ceil(dimCols(2)/2)+1:center(3)+floor(dimCols(2)/2), :) = corY;
panel(dimRows(1)+1:dimRows(1)+dimRows(2), center(4)-ceil(dimCols(2)/2)+1:center(4)+floor(dimCols(2)/2), :) = corZ;
%
panel(dimRows(1)+dimRows(2)+1:dimRows(1)+dimRows(2)+dimRows(3), center(1)-ceil(dimCols(3)/2)+1:center(1)+floor(dimCols(3)/2), :) = traB;
panel(dimRows(1)+dimRows(2)+1:dimRows(1)+dimRows(2)+dimRows(3), center(2)-ceil(dimCols(3)/2)+1:center(2)+floor(dimCols(3)/2), :) = traX;
panel(dimRows(1)+dimRows(2)+1:dimRows(1)+dimRows(2)+dimRows(3), center(3)-ceil(dimCols(3)/2)+1:center(3)+floor(dimCols(3)/2), :) = traY;
panel(dimRows(1)+dimRows(2)+1:dimRows(1)+dimRows(2)+dimRows(3), center(4)-ceil(dimCols(3)/2)+1:center(4)+floor(dimCols(3)/2), :) = traZ;
% Create the gifs
if ~exist([saveDir,'\Panels\'],'dir')
    mkdir([saveDir,'\Panels\']);
end
create_GIF(panel,[saveDir,'\Panels\',saveName,'_panel'], fps, 1);
end



