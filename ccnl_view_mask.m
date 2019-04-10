function ccnl_view_mask(mask)

    % View a mask(s) passed as a 3D array, or cell array of 3D arrays
    % instead of 3D array. Could pass .nii filename(s) instead.
    %
    % Requires MATLAB >= 2018
    %
    % USAGE:
    %   ccnl_check_mask(masks)
    %
    % EXAMPLES:
    %   ccnl_check_mask(ccnl_load_mask('masks/hippocampus.nii'))
    %   ccnl_check_mask('masks/hippocampus.nii')
    %   ccnl_check_mask({'masks/hippocampus.nii', 'masks/mean.nii'})
    %
    % Momchil Tomov, Nov 2016

    if ~iscell(mask)
        mask = {mask};
    end

    masks = mask; 
    filenames = {};
    for i = 1:length(masks)
        if ischar(masks{i})
            filenames = [filenames, masks(i)];
        else
            mask = double(masks{i});

            if ~exist('temp', 'dir')
                mkdir('temp');
            end
            filename = fullfile('temp', sprintf('tmp_%d.nii', i));
            
            filename
            niftiwrite(mask, filename);
            filenames = [filenames, {filename}];
        end
    end
    spm_check_registration(char(filenames));
