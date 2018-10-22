% helper function that extracts quantities from .nii files
%
% Momchil Tomov, Aug 2018
%
function y = get_beta_or_tmap_helper(regressor, modeldir, Vmask, mask, names, mask_format, prefix)
    n = 0;
    for i = 1:length(names)
        if (ischar(regressor) && (~isempty(strfind(names{i},[regressor,'*'])) || ~isempty(strfind(names{i},[regressor,'^'])) || endsWith(names{i},regressor))) ...
           || (~ischar(regressor) && i == regressor)
            n = n + 1;
            V = spm_vol(fullfile(modeldir,sprintf('%s_%04d.nii',prefix,i)));    % residual variance or t-statistic image
            Y = spm_read_vols(V);

            switch mask_format 
                
                case 'inds_or_binary' % a list of voxel indices or a binary vector in NATIVE space
                    y(n,:) = Y(mask);

                case 'cor' % a list of coordinates
                    cor = mask;
                    inds = sub2ind(size(Y),cor(:,1),cor(:,2),cor(:,3));
                    y(n,:) = Y(inds);

                case 'mask' % a 3D mask in MNI or native space
                    if all(size(Y) == size(mask))
                        y(n,:) = Y(mask);
                    else
                        [mask_cor1, mask_cor2, mask_cor3] = ind2sub(size(mask),find(mask==1));
                        cor = mni2cor(cor2mni([mask_cor1 mask_cor2 mask_cor3], Vmask.mat),V.mat);
                        inds = sub2ind(size(Y),cor(:,1),cor(:,2),cor(:,3));
                        y(n,:) = Y(inds);
                    end

                otherwise
                    assert(false, 'Error: incorrect mask or voxel list');
            end
        end
    end
end

