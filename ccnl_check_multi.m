function ccnl_check_multi(EXPT, glmodels, subjs, runs, do_plot)

% Sanity test your multi structure. Useful to run before running ccnl_fmri_glm.
% Makes sure EXPT.create_multi doesn't crash and that the regressors
% are linearly independent.
%
% USAGE:
%   ccnl_check_multi(EXPT, glmodels)
%   ccnl_check_multi(EXPT, glmodels, subjs)
%   ccnl_check_multi(EXPT, glmodels, subjs, runs)
%   ccnl_check_multi(EXPT, glmodels, subjs, runs, do_plot)
%
% INPUT:
%   EXPT = experiment structure
%   glmodels = which models to test, e.g. 1:20
%   subjs (optional) = subject indices to test, e.g. 1:31. Defaults to all subjects
%   runs (optional) = runs to test, e.g. 1:8. Defaults to all runs.
%   do_plot (optional) = whether to plot the regressors convolved with the HRF
%
% EXAMPLES:
%   ccnl_check_multi(exploration_expt(), 1)
%   ccnl_check_multi(exploration_expt(), 1:3)
%   ccnl_check_multi(exploration_expt(), 1, 1:31)
%   ccnl_check_multi(exploration_expt(), 1, 1:31, 1)
%   ccnl_check_multi(exploration_expt(), 1, 1, 1:8)
%   ccnl_check_multi(exploration_expt(), 1, 1, 1, true)
% 

% set default parameters
%
if nargin < 5
    do_plot = false;
end

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
            [X, names] = ccnl_get_design(EXPT, glmodel, subj, run);

            % optionally plot regressors
            if do_plot
                for i = 1:size(X,2)
                    subplot(size(X,2), 1, i);
                    plot(X(:,i));
                    legend(names(i), 'Interpreter', 'none');
                end
            end

            % check for dependent columns
            [R, jb] = rref(X);
            if numel(jb) < size(X,2)
                disp(['Found co-linear regressors for glmodel ', num2str(glmodel), ', subj ', num2str(subj), ', run ', num2str(run), ':']);

                % print dependent columns
                c = logical(ones(1,size(X,2)));
                c(jb) = 0;
                c = find(c);
                for i = c
                    s = [names{i}, ' = '];
                    first = true;
                    for j = 1:size(R,1)
                        if abs(R(j,i)) > 1e-15
                            if ~first
                                s = [s, ' + '];
                            end
                            s = [s, num2str(R(j,i)), ' * ', names{j}];
                            first = false;
                        end
                    end
                    disp(s);
                end
            end
        end
    end
end
