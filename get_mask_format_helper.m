% helper function that figures out the what format a lit of voxels is in
% potentially actually loads the mask
%
% Momchil Tomov, Aug 2018
%
function [mask_format, mask, Vmask] = get_mask_format_helper(mask)

    % figure out how voxels are provided
    mask_format = 'none';
    if ischar(mask)
        % load mask
        Vmask = spm_vol(mask);
        Ymask = spm_read_vols(Vmask);
        mask = Ymask~=0 & ~isnan(Ymask);
        mask_format = 'mask';
    else
        Vmask = [];
        if ndims(mask) == 3
            mask_format = 'mask';
        elseif ndims(mask) == 2
            if size(mask,2) == 3
                mask_format = 'cor';
            elseif size(mask,1) == 1 || size(mask,2) == 1
                mask_format = 'inds_or_binary';
            end
        end
    end
