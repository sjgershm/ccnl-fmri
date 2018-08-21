function ccnl_plot_regressors(EXPT, glmodel, subj, run)

% Plot the regressors for a given glmodel, subject & run.
% Requires the SPM structure to have been generated i.e. the single-subject
% GLM should have been run (ccnl_fmri_glm).
%
% USAGE:
% ccnl_plot_regressors(EXPT, glmodel, subj, run)
%
% INPUTS:
%   EXPT - experiment structure
%   glmodel - model number
%   subj - subject id
%   run = which run
%
% EXAMPLE:
% ccnl_plot_regressors(exploration_expt(), 1, 31, 8)
%
% Momchil Tomov, Aug 2018

include_motion = false;
%mask = 'hippocampus.nii';

multi = EXPT.create_multi(glmodel, subj, run);

modeldir = fullfile(EXPT.modeldir,['model',num2str(glmodel)],['subj',num2str(subj)]);
load(fullfile(modeldir,'SPM.mat'));

TR = EXPT.TR;
sess_prefix = ['Sn(', num2str(run), ')'];
trs = 1 : TR : TR*length(SPM.Sess(run).row); % or start at 0? how does slice timing interpolation work in SPM?


figure;

% which regressors to display
%
cols = SPM.Sess(run).col;
if ~include_motion
    cols = cols(1:end-6); % ditch motion regressors
end


% iterate over regressors
%
plot_idx = 0;
for i = cols
    assert(strncmp(SPM.xX.name{i}, sess_prefix, length(sess_prefix)));
    
    plot_idx = plot_idx + 1;
    subplot(length(cols) + 1, 1, plot_idx);
    hold on;
    h = [];
    
    % plot trial onsets / offsets as vertical dashed lines
    %
    %feedback_onsets = cellfun(@str2num, data.actualChoiceOnset(which_rows)');
    %for t=feedback_onsets
    %    plot([t t], [-1 1], '--', 'Color', [0.8 0.8 0.8]);
    %end
    
    % plot original regressor from model
    %
    eps = 1e-6;
    leg = {}; % legend
    for j = 1:length(multi.names)
        n = length(multi.onsets{j});
        onsets = multi.onsets{j};
        durations = multi.durations{j};
        if ~(size(durations, 1) == size(onsets, 1) && size(durations, 2) == size(onsets, 2))
            % need to rotate one of the vectors
            durations = durations';
        end
        assert(size(durations, 1) == size(onsets, 1) && size(durations, 2) == size(onsets, 2));
        x = [onsets - eps; onsets; onsets + durations; onsets + durations + eps];
        x = x(:)';
        x = [0 x max(trs)];
        if ~isempty(strfind(SPM.xX.name{i}, [' ', multi.names{j}, '*']))
            y = [zeros(1,n); ones(1,n); ones(1,n); zeros(1,n)];
            y = y(:)';
            y = [0 y 0];
            h = [h, plot(x, y, 'LineWidth', 1, 'Color', 'red')];
            leg = [leg; multi.names{j}];
        end
            
        if isfield(multi, 'pmod') && j <= length(multi.pmod)
            for k = 1:length(multi.pmod(j).name)
                if ~isempty(strfind(SPM.xX.name{i}, ['x', multi.pmod(j).name{k}, '^']))
                    y = reshape(multi.pmod(j).param{k}, 1, n);
                    y = [zeros(1,n); y; y; zeros(1,n)];
                    y = y(:)';
                    y = [0 y 0];
                    h = [h, plot(x, y, 'LineWidth', 1, 'Color', 'red')];
                    leg = [leg; ['pmod: ', multi.pmod(j).name{k}]];
                end
            end
        end
    end

    
    % plot regressor convolved with HRF
    %
    h = [h, plot(trs, SPM.xX.X(SPM.Sess(run).row, i)', 'Color', 'blue')];
    leg = [leg; {SPM.xX.name{i}}];
    
    % TODO plot beta
    %
    %beta_vec = ccnl_get_beta(EXPT, glmodel, i, mask, [subj]);
    
    yL = get(gca,'YLim');
    ylim([yL(1), yL(2) + 0.1]);    
    title(SPM.xX.name{i}, 'Interpreter', 'none');
    legend(h, leg, 'Interpreter', 'none');
        
    hold off
end
