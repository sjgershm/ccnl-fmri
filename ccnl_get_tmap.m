function tmap = ccnl_get_tmap(EXPT,model,regressor,mask,subjects)
    
    % Extract t-statistics from a mask. Also see ccnl_get_beta.
    %
    % USAGE: tmap = ccnl_get_tmap(EXPT,model,regressor,mask,[subjects])
    %
    % INPUTS:
    %   EXPT - experiment structure
    %   model - model number
    %   regressor - regressor name or regressor index
    %   mask - a mask image name (e.g., 'mask.nii'), a set of voxel
    %          indices, or a binary vector
    %   subjects (optional) - which subjects to analyze (default all subjects)
    %
    % OUTPUTS:
    %   tmap - [nSubjects x nVoxels] t-statistics
    %
    % Momchil Tomov, Jul 2017
    
    if nargin < 4; subjects = 1:length(EXPT.subject); end
    
    % load mask
    if ischar(mask)
        Vmask = spm_vol(mask);
        Ymask = spm_read_vols(Vmask);
        mask = Ymask~=0 & ~isnan(Ymask);
    end

    if ischar(regressor)
        regressor = strtrim(regressor);
    end

    % load contrasts
    load(fullfile(EXPT.modeldir,['model',num2str(model)],'contrasts.mat'));

    for s = 1:length(subjects)
        subj = subjects(s);
        modeldir = fullfile(EXPT.modeldir,['model',num2str(model)],['subj',num2str(subj)]);

        t = get_beta_or_tmap_helper(regressor, modeldir, mask, contrasts, 'spmT');
    
        if size(t, 1) > 1
            load(fullfile(modeldir,'SPM.mat'));
            b = get_beta_or_tmap_helper(regressor, modeldir, mask, SPM.xX.name, 'beta');
            assert(isequal(size(b), size(t)), 'Should have 1 t-map only or exactly 1 t-map per regressor');

            se = b ./ t;
            se = sqrt(nansum(se.^2)) / size(se, 1);
            b = nanmean(b, 1);
            t = b ./ se;
        end
        tmap(s,:) = t;
        
    end

end


% helper function
%
function y = get_beta_or_tmap_helper(regressor, modeldir, mask, names, prefix)
    n = 0;
    for i = 1:length(names)
        if ischar(regressor) && (~isempty(strfind(names{i},[regressor,'*'])) || ~isempty(strfind(names{i},[regressor,'^'])) || endsWith(names{i},regressor)) ...
           || ~ischar(regressor) && i == regressor
            n = n + 1;
            V = spm_vol(fullfile(modeldir,sprintf('%s_%04d.nii',prefix,i)));    % residual variance or t-statistic image
            Y = spm_read_vols(V);
            
            if all(size(Y) == size(mask))
                y(n,:) = Y(mask);
            else
                [mask_cor1, mask_cor2, mask_cor3] = ind2sub(size(mask),find(mask==1));
                cor = mni2cor(cor2mni([mask_cor1 mask_cor2 mask_cor3], Vmask.mat),V.mat);
                inds = sub2ind(size(Y),cor(:,1),cor(:,2),cor(:,3));
                y(n,:) = Y(inds);
            end
        end
    end
end

