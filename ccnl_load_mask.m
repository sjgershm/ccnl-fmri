function [mask, V, Y] = ccnl_load_mask(mask_filename)

    % Load mask from .nii file as a 3D array.
    %
    % INPUT:
    %   mask_filename = path to .nii file with mask
    %
    % OUTPUT:
    %   mask = 3D binary mask for the non-0 non-nan voxels from Ymask
    %   V = volume struct from spm_vol, could be used to save the mask with spm_write_vol(V, Y) or spm_write_vol(V, mask). IMPORTANT: make sure to change V.fname first! Otherwise, will overwrite original file
    %   Y = mask loaded with from spm_read_vols 
    %
    % USAGE:
    %   [mask, V, Y] = ccnl_load_mask(mask_filename)
    %
    % EXAMPLES:
    %   [mask, V, Y] = ccnl_load_mask('masks/hippocampus.nii')
    %
    V = spm_vol(mask_filename);
    Y = spm_read_vols(V);
    mask = Y~=0 & ~isnan(Y);
