function [X, names, X_raw, res] = ccnl_get_design(EXPT, glmodel, subj, run)

    % Compute the design matrix for a given run by convolving the 
    % regressors with the HRF. Useful for plotting and sanity checks.
    % Notice that this does not rely on the SPM structure (unlike ccnl_plot_regressors).
    %
    % USAGE:
    %   ccnl_get_design(EXPT, glmodel, subj, run)
    %
    % INPUT:
    %   EXPT = experiment structure
    %   glmodels = which model
    %   subj = which subject id
    %   run = which run
    %
    % OUTPUT:
    %   X = [nSeconds x nRegressors] design matrix for given run (assuming TR = 1 s)
    %   names = [1 x nRegressors] cell array of regressor names
    %   X_raw = [nSeconds*res x nRegressors] raw regressors, not convolved & subsampled
    %   res = resolution (ms)
    %
    % EXAMPLE:
    %   [X, names] = ccnl_get_design(exploration_expt(), 18, 1, 1)
    % 
    % Momchil Tomov, Aug 2018


    res = 1000; % ms res
    multi = EXPT.create_multi(glmodel, subj, run);

    % find max onset
    run_dur = EXPT.run_duration * res; % s => ms

    % iterate over events
    names = {};
    for j = 1:length(multi.names)
        onsets = multi.onsets{j} * res;
        durations = multi.durations{j} * res;

        x = populate_regressor(onsets, durations, ones(size(onsets)), run_dur);
        if j == 1
            X_raw = x;
        else
            X_raw = [X_raw x];
        end
        x = convolve_and_subsample(x, res);
        if j == 1
            X = x;
        else
            X = [X x];
        end
        names = [names, multi.names{j}];

        % iterate over pmods for event
        if isfield(multi, 'pmod') && j <= length(multi.pmod)
            for k = 1:length(multi.pmod(j).name)
                assert(length(onsets) == length(multi.pmod(j).param{k}), ['multi.pmod(', num2str(j), ').param{', num2str(k), '} has the wrong number of elements']);

                x = populate_regressor(onsets, durations, multi.pmod(j).param{k}, run_dur);
                X_raw = [X_raw x];
                x = convolve_and_subsample(x, res);
                X = [X x];
                names = [names, multi.pmod(j).name(k)];
            end
        end
    end

end


function x = populate_regressor(onsets, durations, values, max_time)
    x = zeros(ceil(max_time),1);
    for i = 1:length(onsets)
        x(round(onsets(i)):round(min(max_time, onsets(i)+durations(i)))) = values(i);
    end
end



function x = convolve_and_subsample(x, res)
    hrf = spm_hrf(1 / res);
    x = conv(x, hrf);
    x = x(1:res:end);
end
