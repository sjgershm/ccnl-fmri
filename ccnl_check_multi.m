function ccnl_check_multi(EXPT, glmodels, subjs, runs)

% Sanity test your multi structure. Useful to run before running ccnl_fmri_glm.
%
% USAGE:
%   ccnl_check_multi(EXPT, glmodels)
%   ccnl_check_multi(EXPT, glmodels, subjs)
%   ccnl_check_multi(EXPT, glmodels, subjs, runs)
%
% INPUT:
%   EXPT = experiment structure
%   glmodels = which models to test, e.g. 1:20
%   subjs (optional) = subject indices to test, e.g. 1:31. Defaults to all subjects
%   runs (optional) = runs to test, e.g. 1:8. Defaults to all runs.
%
% EXAMPLES:
%   ccnl_check_multi(exploration_expt(), 1)
%   ccnl_check_multi(exploration_expt(), 1:3)
%   ccnl_check_multi(exploration_expt(), 1, 1:31)
%   ccnl_check_multi(exploration_expt(), 1, 1, 1)
%   ccnl_check_multi(exploration_expt(), 1, 1, 1:8)
% 

% set default parameters
%
if nargin < 4 
    all_runs = true;
else
    all_runs = false;
end

if nargin < 3
    subjs = 1:length(EXPT.subject);
end

for glmodel = glmodels
    for subj = subjs
        if all_runs
            runs = 1:length(EXPT.subject(subj).functional);
        end
        for run = runs
            multi = EXPT.create_multi(glmodel, subj, run);
        end
    end
end
