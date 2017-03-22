function beta = ccnl_get_beta(EXPT,model,regressor,mask,subjects)
    
    % Extract beta coefficients from a mask.
    %
    % USAGE: beta = ccnl_get_beta(EXPT,model,regressor,mask,[subjects])
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
    %   beta - [nSubjects x nVoxels] beta coefficients
    %
    % Sam Gershman, Nov 2016
    
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

    for s = 1:length(subjects)
        subj = subjects(s);
        modeldir = fullfile(EXPT.modeldir,['model',num2str(model)],['subj',num2str(subj)]);
        load(fullfile(modeldir,'SPM.mat'));
       
        n = 0; clear y
        for i = 1:length(SPM.xX.name)
            if ischar(regressor) && (~isempty(strfind(SPM.xX.name{i},[regressor,'*'])) || ~isempty(strfind(SPM.xX.name{i},[regressor,'^']))) ...
               || ~ischar(regressor) && i == regressor
                n = n + 1;
                V = spm_vol(fullfile(modeldir,sprintf('beta_%04d.nii',i)));    % residual variance image
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
        
        beta(s,:) = nanmean(y, 1);
        
    end
