function [vifs, names] = ccnl_vifs(EXPT, model)

    % Compute Variance Inflation Factors (VIFs) for given GLM.
    % Need to run ccnl_fmri_glm first to compute the SPM.mat's
    %
    % USAGE:
    %   [vifs, names] = ccnl_vifs(EXPT, glmodel)
    %
    % INPUTS:
    %   EXPT - experiment structure
    %   model - model number
    %
    % OUTPUTS:
    %   vifs - [nRegressors] variance inflation factors
    %   names - [nRegressors] regressor names
    %

    subjects = 1:length(EXPT.subject);

    for s = 1:length(subjects)
        subj = subjects(s);
        modeldir = fullfile(EXPT.modeldir,['model',num2str(model)],['subj',num2str(subj)]);
        if ~exist(fullfile(modeldir,'SPM.mat'), 'file')
            % keep looking; we only need to find one subject with a SPM.mat
            continue;
        end

        out = scn_spm_design_check(modeldir, 'events_only');

        break; % we only need to find one subject with a SPM.mat
    end

    if ~exist('out', 'var')
        assert(false, 'No SPM.mat found for any subject -- make sure to run GLM first.');
    end
    names = out.name;
    vifs = out.allvifs;

end




function out = scn_spm_design_check(spm_results_dir, varargin)

% shamelessly stolen from https://github.com/canlab/CanlabCore/
%
% Run in a single-subject (first-level) SPM directory to check 
% design matrix variance inflation factors and high-pass filtering.
% Prints out table of regressors and their above-threshold VIFs (see options).
% Saves .png images of the key figures.
%
% :Usage:
% ::
%
%     scn_spm_design_check(spm_results_dir, varargin)
%
% :Optional Inputs:
%
%   **'events_only':**
%        Show plots and diagnostics for ONLY events, not nuisance covariates or
%        other user-specified regressors.  Useful when you have many nuisance
%        covs.
%
%   **'vif_thresh', t':**
%        Only regressors with a VIF > t will be printed in VIF table.
%
%   **'sort_by_vif'':**
%        Sort regressors in VIF table by VIF (DEFAULT: order regressors as in model).   
%
% Calls: scn_spm_choose_hpfilter.m, scn_spm_get_events_of_interest.m
%
% :Examples:
% ::
%
%    scn_spm_design_check(pwd, 'events_only');
%
% ..
%    Updated: Tor Wager, Aug 2010; Oct 2011: Add 'events_only'; July 2012:
%    fixed for parametric modulators. Luka Ruzic, Sept 2012: added VIF tables.
%    Wani Woo, Apr, 2018: added an output (out) to return vif values 
% ..


if nargin < 1, spm_results_dir = pwd; end

spmfilename = fullfile(spm_results_dir, 'SPM.mat');
if ~exist(spmfilename, 'file')
    error('SPM.mat does not exist in %s\n', spm_results_dir); 
end
load(spmfilename);

VIFTHRESH = 1.3;
EVENTS_ONLY = false;
SORTBYVIF = false;

i=1;
while i<=numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'events_only'
                EVENTS_ONLY = true;
            case 'vif_thresh'
                i=i+1;
                VIFTHRESH = varargin{i};
            case 'sort_by_vif'
                SORTBYVIF = true;
            case 'sort_by_reg'
                % backwards compatibility; do nothing
            otherwise
                error(['Unrecognized argument: ' varargin{i}])
        end
    else
        error(['Unrecognized argument: ' varargin{i}])
    end
    i=i+1;
end


%% Optional inputs


% Gets events of interest: All regressors, or events only if 'events_only'
% is input as keyword
if EVENTS_ONLY
    wh_cols = scn_spm_get_events_of_interest(SPM, 'events_only');
else
    wh_cols = scn_spm_get_events_of_interest(SPM);
end


%%

%create_figure('X', 2, 2); 
%imagesc(zscore(SPM.xX.X)); set(gca, 'YDir', 'Reverse');
%colormap gray
%title('Full Design Matrix (zscored)');
%axis tight
%drawnow
%
%subplot(2, 2, 2)
%imagesc(zscore(SPM.xX.X(:, wh_cols))); set(gca, 'YDir', 'Reverse');
%colormap gray
%title('Design Matrix: Of-interest (zscored)');
%axis tight
%drawnow

% Variance Inflation Factors for regs of interest (iC)

warning off % some nuisance covs could be redundant; we don't care.
allvifs = getvif(SPM.xX.X(:, SPM.xX.iC), 0);
warning on

%subplot(2, 2, 3);
allvifs = allvifs(wh_cols);

%plot(allvifs, 'ko', 'MarkerFaceColor', [1 .5 0]);
%
%ylabel('Variance inflation factor'); xlabel('Predictor number');
%plot_horizontal_line(1, 'k');
%plot_horizontal_line(2, 'b--');
%plot_horizontal_line(4, 'r--');
%plot_horizontal_line(8, 'r-');
%disp('Variance inflation: 1 (black line) = minimum possible (best)');
%disp('Successive lines indicate doublings of variance inflation factor.');
%title('Var. Inflation (VIFs) in full design');
%
%% Variance Inflation Factors for ONLY of-interest
%% now we care if nuisance covs are redundant.
%vifs = getvif(SPM.xX.X(:, wh_cols), 0);
%subplot(2, 2, 4);
%plot(vifs, 'ko', 'MarkerFaceColor', [1 .5 0]);
%
%ylabel('Variance inflation factor'); xlabel('Predictor number');
%plot_horizontal_line(1, 'k');
%plot_horizontal_line(2, 'b--');
%plot_horizontal_line(4, 'r--');
%plot_horizontal_line(8, 'r-');
%title('VIFs for ONLY of-interest regs');

% table
if SORTBYVIF
    [ignore ord] = sort(allvifs,'descend'); %#ok
else
    ord = 1:numel(allvifs);
end
fprintf('\nVIFs greater than %d\n',VIFTHRESH)
fprintf('\n%7s   %s\n','VIF','REGRESSOR')
for i=1:numel(ord)
    if allvifs(ord(i)) >= VIFTHRESH
        fprintf('%7.2f   reg%04d: "%s"\n',allvifs(ord(i)),wh_cols(ord(i)),SPM.xX.name{wh_cols(ord(i))})
    end
end

%try 
%    scn_export_papersetup; 
%catch
%end
%saveas(gcf, 'Variance_Inflation', 'png');
%disp('Saved Variance_Inflation.png in SPM directory');

%if EVENTS_ONLY
%    scn_spm_choose_hpfilter(spm_results_dir, 'events_only');
%else
%    scn_spm_choose_hpfilter(spm_results_dir);
%end

%spm_efficiency('SPM.mat');
% saveas(gcf, 'SPM_efficiency', 'png');
% disp('Saved SPM_efficiency.png in current directory');

%try 
%    scn_export_papersetup(500); 
%catch
%end
%saveas(gcf, 'High_pass_filter_analysis', 'png');
%disp('Saved High_pass_filter_analysis.png in SPM directory');

out.allvifs = allvifs;
out.name = SPM.xX.name(wh_cols);

end





function vif = getvif(X, no_add_intercept, varargin)

% function vif = getvif(model design matrix (X), [no_add_intercept], varargin)
%
% Optional arguments:
%
% 1.  first varargin must be flag for no_add_intercept.
% The only case where we may not want to add an intercept is if we already
% have intercepts for each subset of observations (e.g., each run in fmri)
% and these are sufficient.
%
% 2.
% pass in 'wh' followed by vector of indices to return the vifs for only those columns,
% though vif will still be calculated on entire model.
%
% refactored and updated documentation: by tor wager, aug 2015
%
% Examples:
% Generate X:
% [X, delta, delta_hires, hrf] = onsets2fmridesign(ons, TR, scanLength, 'hrf', 'parametric_standard', pm);
% OR
% X = mvnrnd([1 1 1 1], eye(4), 100);
%
% vifs = getvif(X)
%
% getvif(X, 0, 'plot')
%
%
% See also:
% scn_spm_design_check

% Remove intercept, we will add later to each model
X = intercept(X, 'remove');         % remove if present

[n, k] = size(X);

% Optional arguments
% -------------------------------------------------------------------------

if nargin < 2, no_add_intercept = 0; end

ind = strcmp('wh', varargin);

if ~isempty(ind) && sum(ind) ~= 0
    
    wh=varargin{ind+1};
    if size(wh, 1) > size(wh, 2), wh=wh';end %transpose if needed
    
else
    
    wh = 1:k;
    
end

doplot = 0;

if any(strcmp('plot', varargin))
    doplot = 1;
end

% calculate
% -------------------------------------------------------------------------

% initialize, exclude the intercept
vif = [];  %ones(1, length(wh));

intcpt = ones(n, 1);

for i = wh
    
    if no_add_intercept
        Xi = X;
    else
        Xi = [X intcpt];
    end
    
    % check that model has an intercept, or that all preds are
    % mean-centered.  if both of these conditions are false, Rsquared is uninterpretable see [http://statistiksoftware.blogspot.com/2013/01/why-we-need-intercept.html]
    % and thus VIF is uninterpretable.  Yoni, May 2015
    hasint = any(var(Xi)==0);
    totallymc = sum(mean(Xi))==0;
    
    if ~hasint && ~totallymc
        warning('Model has no intercept, and model is not mean-centered; VIF is not interpertable');
        vif = NaN * zeros(1, length(wh));
        return
    end
    
    
    y = Xi(:, i);
    
    Xi(:, i) = [];
    
    b = Xi\y;
    fits = Xi * b;
    
    rsquare = var(fits) / var(y);
    
    % sometimes rounding error makes rsquare>1
    if rsquare >= 1, rsquare = .9999999; end
    
    vif(end + 1) = 1 / (1 - rsquare);
    
    
end  % regressor

if doplot
    
    hold on;
    plot(wh, vif, 'ko', 'MarkerFaceColor', [1 .5 0], 'MarkerSize', 8);
    
    ylabel('Variance inflation factor'); xlabel('Predictor number');
    plot_horizontal_line(1, 'k');
    h = plot_horizontal_line(2, '--'); set(h, 'Color',[0 .3 .7]);
    h = plot_horizontal_line(4, '--'); set(h, 'Color', [0 .7 .3]);
    h = plot_horizontal_line(8, '-'); set(h, 'Color', [.7 .3 0]);
    
    disp('Variance inflation: 1 (black line) = minimum possible (best)');
    disp('Successive lines indicate doublings of variance inflation factor.');
    disp('Red boxes have extremely high VIFs, perfect multicolinearity');
    
    title('Var. Inflation (VIFs) in full design');
    
    mymax = max(10, max(vif(vif < 1000))+3);
    
    set(gca, 'YLim', [0 mymax], 'FontSize', 20, 'XTick', 1:length(wh), 'XLim', [min(wh)-.5 max(wh) + .5]);
    
    whhigh = find(vif > 1000);
    for i = 1:length(whhigh)
        
        drawbox(wh(whhigh(i)) - .3, .6, 0, mymax, [1 0 0]);
        
    end
    
end

end % function



function X = intercept(X, meth)
% Intercept-related functions for working with design matrices, etc.
%
% Return which columns are intercept(s), if any
% ::
%
%    wh = intercept(X, 'which')
%
% Remove an intercept, if there is one
% ::
%
%    X = intercept(X, 'remove');
%
% Add an intercept to the end
% ::
%
%    X = intercept(X, 'add');
%
% Ensure that the intercept is at the end, moving or adding it as necessary
% ::
%
%    X = intercept(X, 'end');
%

wh_is_intercept = find( ~any(diff(X)) ); % which column is the intercept?

switch meth
    
    case 'which'
        X = wh_is_intercept;
        
    case 'remove'
        
        if ~isempty( wh_is_intercept ), X(:, wh_is_intercept) = []; end
        
    case 'add'
        
        X(:, end+1) = 1;
        
    case 'end'
        
        X = intercept(X, 'remove');
        X = intercept(X, 'add');
        
    otherwise
        
        error('Unknown method');
        
end

end



function wh_cols = scn_spm_get_events_of_interest(SPM, varargin)
% Gets events of interest.
%
% :Usage:
% ::
%
%     wh_cols = scn_spm_get_events_of_interest(SPM, varargin)
%
% All regressors, or events only if 'events_only' is input as keyword
% 'from_multireg':  followed by an integer, to include first n columns from
% the multireg R matrix as "of interest".  only works with 'events_only'
% flag, of course.

events_only = false;  % false for all regressors, true for events only
numcols_to_add_from_multireg = 0;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % reserved keywords
            case 'events_only', events_only = true;
            case 'from_multireg', numcols_to_add_from_multireg = varargin{i+1};
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end


if events_only
    
    nsess = length(SPM.Sess);
    
    % Tor's code to get only events
    for i = 1:nsess
        ncols_in_run(i, 1) = length(SPM.Sess(i).col);
        
        %n_ofinterest(i) = length(SPM.Sess(i).U);
        % Modified to handle parametric modulators as well. 7/2012
        n_ofinterest(i) = length(cat(2, SPM.Sess(i).U.name));
        
    end
    
    if ~all(n_ofinterest == n_ofinterest(1))
        disp('Warning! Different numbers of events of interest in each run.');
        disp('If this is what you intended, OK, but will use max val and so will include some of no interest')
    end
        
    firstcol = cumsum([1; ncols_in_run(1:end-1)]);
    
    for i = 1:max(n_ofinterest) + numcols_to_add_from_multireg
        wh_cols(:, i) = firstcol + (i - 1);
    end
    
    wh_cols = wh_cols';
    wh_cols = wh_cols(:);
    
    wh_cols = SPM.xX.iC(wh_cols);  % iC should be index of all...so redundant, but include for logical consistency
else
    wh_cols = SPM.xX.iC;
end

end % function

%%





