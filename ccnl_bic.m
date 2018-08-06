function bic = ccnl_bic(EXPT,model,mask,subjects)

% Compute Bayesian information criterion (BIC) for a model.
%
% USAGE: bic = ccnl_bic(EXPT,model,mask,[subjects])
%
% INPUTS:
%   EXPT - experiment structure
%   model - model number
%   mask - a mask image name (e.g., 'mask.nii'), a set of voxel
%          indices, or a binary vector
%   subjects (optional) - which subjects to analyze (default all subjects)
%
% OUTPUTS:
%   bic - [nSubjects x 1] vector of BIC values
%
% Sam Gershman, Oct 2016

if nargin < 3; subjects = 1:length(EXPT.subject); end

% load mask
if ischar(mask)
    Vmask = spm_vol(mask);
    Ymask = spm_read_vols(Vmask);
    mask = Ymask~=0 & ~isnan(Ymask);
end

bic = zeros(length(subjects),1);
for s = 1:length(subjects)
    subj = subjects(s);
    modeldir = fullfile(EXPT.modeldir,['model',num2str(model)],['subj',num2str(subj)]);
    load(fullfile(modeldir,'SPM.mat'));
    [N,K] = size(SPM.xX.X); % number of datapoints and parameters
    V = spm_vol(fullfile(modeldir,'ResMS.nii'));    % residual variance image
    Y = spm_read_vols(V);
    
    if all(size(Y) == size(mask))
        bic(s) = N*nansum(log(Y(mask))) + K*log(N);
       
    else
        [mask_cor1, mask_cor2, mask_cor3] = ind2sub(size(mask),find(mask==1));
        res_cor = mni2cor(cor2mni([mask_cor1 mask_cor2 mask_cor3], Vmask.mat),V.mat);
        res_inds = sub2ind(size(Y),res_cor(:,1),res_cor(:,2),res_cor(:,3));
        disp(sum(isnan(Y(res_inds))))
        bic(s) = N*nansum(log(Y(res_inds))) + K*log(N);
    end
end
