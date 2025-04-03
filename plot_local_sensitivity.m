% postprocess output from all the local sensitivity analysis

incr = 1; % increase (1) or decrease (0) parameters by 5%

if incr
    type_change = 'I';
else
    type_change = 'D';
end
fname = strcat('./Rat_Data/localsensitivity_male_',type_change,'.mat');
male_dat = load(fname);

[male_xlabs, male_vals] = getdata(male_dat);

xlabels = male_xlabs;

% Make figures
temp = male_vals;

ylabels = categorical(["RSNA"; "GFR"; "MAP"; "[ALD]"; "PRA"; "[PTH]"; ...
                       "[1,25(OH)_2D_3]"; "[Mg^{2+}]"; "[Ca^{2+}]"]);

fsize = 14;
colmap = parula;
cmissdat = 'w';
labmissdat = '<1.0%';

% remove Nan values
[male_vals_rm, xlabs_rm, ylabs_rm] = remove_all_nan_rows_and_cols(male_vals, xlabels, ylabels);

% figure with removing 
fig_remove = 1;
if fig_remove
f2 = figure('Position', [100, 100, 1000, 1500]);
f2.Position = [1 1 3000 600];
clf

h1 = heatmap(xlabs_rm, ylabs_rm, male_vals_rm,...
                'colormap', colmap,...
                'MissingDataColor', cmissdat, 'MissingDataLabel', labmissdat, ...
                'colorbarvisible', 'on');
h1.FontSize = fsize;
s = struct(h1);
s.XAxis.TickLabelRotation = 90;
s.XAxis.FontSize = 16;
s.YAxis.FontSize = 16;
end

%------------------
% functions used
%------------------
function [xlabels, round_data] = getdata(dat)
    frac_sens = dat.frac_change;
    xlabels = dat.pars_name;
    round_data = round(frac_sens', 2, 'significant');
    [r,c] = find(abs(round_data) < 1.0);% r - row values, c - column value
    for ii = 1:length(r)
        round_data(r(ii),c(ii)) = NaN;
    end    
end

function [male_vals_rm, xlabs_rm, ylabs_rm] = remove_all_nan_rows_and_cols(male_vals, xlabels, ylabels)
    % Find columns with all NaN values
    all_nan_cols = all(isnan(male_vals));
    kept_col_indices = find(all_nan_cols == 0);

    % Remove rows and columns with all NaN values
    male_vals_rm = male_vals(:, ~all_nan_cols);
    xlabs_rm = xlabels(kept_col_indices);
    ylabs_rm = ylabels;
end