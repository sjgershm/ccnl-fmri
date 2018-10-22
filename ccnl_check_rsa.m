function ccnl_check_rsa(EXPT, rsa_idxs, subjs)

% Sanity test your rsa structure. Useful to run before running ccnl_rsa_searchlight.
% Makes sure EXPT.create_rsa doesn't crash and that the behavioral RDMs can be computed.
%
% USAGE:
%   ccnl_check_rsa(EXPT, rsa_idxs [, subjs])
%
% INPUT:
%   EXPT = experiment structure
%   rsa_idxs = which rsa's to test, e.g. 1:20
%   subjs (optional) = subject indices to test, e.g. 1:31. Defaults to all subjects
%
% EXAMPLES:
%   ccnl_check_rsa(exploration_expt(), 1)
%   ccnl_check_rsa(exploration_expt(), 1:3)
%   ccnl_check_rsa(exploration_expt(), 1, 1:31)
% 
% Momchil Tomov, Oct 2018

if ~exist('subjs', 'var')
    subjs = 1:length(EXPT.subject);
end

% sanity check rsa struct
for rsa_idx = rsa_idxs
    fprintf('rsa_idx = %d\n', rsa_idx);
    for subj = subjs
        fprintf('  subj = %d\n', subj);

        % make sure RSA can be created
        rsa = EXPT.create_rsa(rsa_idx, subj);

        % make sure mask exists
        assert(exist(rsa.mask, 'file') ~= 0, sprintf('Mask %s does not exist', rsa.mask));

        % make sure features and betas match in number
        for i = 1:length(rsa.model)
            assert(sum(rsa.which_betas) == size(rsa.model(i).features, 1), sprintf('Number of betas sum(rsa.which_betas) = %d is different from number of features in rsa.model(%d) == %d', sum(rsa.which_betas), i, size(rsa.model(i).features, 1)));
        end

        % make sure rsa.which_betas is matches GLM betas
        modeldir = fullfile(EXPT.modeldir,['model',num2str(rsa.glmodel)],['subj',num2str(subj)]);
        load(fullfile(modeldir,'SPM.mat'));
        which = contains(SPM.xX.name, rsa.event); % betas for given event
        assert(sum(which) == length(rsa.which_betas), sprintf('There are %d regressors which contain rsa.event ("%s") in their names, but there are %d elements in rsa.which_betas. These should be matching.', sum(which), rsa.event, length(rsa.which_betas)));

        which(which) = rsa.which_betas; % of those, only betas for given trials
    end
end

% try to actually create the behavioral RDMs
for rsa_idx = rsa_idxs
    [Behavioral, control] = ccnl_behavioral_rdms(EXPT, rsa_idx);
end
