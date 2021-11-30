function Vq = hcp_resample_to_mni(vol, x, y, z, gridstep, interp)

conf = mvarconfig('mvarconfig.ini');

if nargin < 6, interp = 'linear'; end
if nargin < 5, gridstep = 8; end

filename = fullfile(conf.DIR_EXTERNAL, ['osl/std_masks/MNI152_T1_' num2str(gridstep) 'mm_brain.nii.gz']);

%[~,dims,scales] = read_avw(filename);
load('/home/chiggins/Documents/MATLAB/MNI52BrainDimsScales');

mni_grid_x =   90 - (0:dims(1)-1)*scales(1);
mni_grid_y = -126 + (0:dims(2)-1)*scales(2);
mni_grid_z =  -72 + (0:dims(3)-1)*scales(3);

[mni_grid_x,mni_grid_y,mni_grid_z,mni_grid_t] = ndgrid(mni_grid_x,mni_grid_y,mni_grid_z, (1:size(vol,4)));

if size(vol,4) > 1,
    Vq = interpn(x,y,z, (1:size(vol,4)),...
                 vol, ...
                 mni_grid_x, ...
                 mni_grid_y, ...
                 mni_grid_z, ...
                 mni_grid_t, interp);
else
    Vq = interpn(x,y,z,...
                 vol, ...
                 mni_grid_x, ...
                 mni_grid_y, ...
                 mni_grid_z, interp);
end%if
