%% simulation input
do_sim_Mg = 0; % change to 1 to perform dietary Mg restriction simulation
do_sim_Ca = 0; % change to 1 to perform dietary Ca restriction simulation
do_sim_vitD3 = 0; % change to 1 to perform 25(OH)D inhibition simulation
do_fig = 1; % change to 1 to plot the results of drug intervention

% simulating 70% inhibition
I_change = 0.3;

tchange=0;
tspan = [0 43800]; % 1 month
[htn_rsna, htn_renin, htn_raa, htn_ald] = deal(1.3,1.3,1.0,1.0);

%% 50% dietary Mg restriction simulation
if do_sim_Mg
    % read baseline datafile
    fname = './Rat_Data/rat_male_ss_data_scenario_normal_combined.mat';
    x1 = load(fname).SSdata;

    % parameters
    species = 'rat';
    sex = 'male';
    [pars_BP, pars_Mg, M] =  get_params_and_mass_matrix(x1,species,sex,htn_rsna,htn_raa,htn_renin,htn_ald);
        
    % options for ode
    options = odeset('Mass',M,'RelTol',1e-2,'AbsTol',1e-3*ones(1,length(x1)));

    % inhibiting dietary Mg intake by 70%
    base_IMg = pars_Mg(28);
    pars_Mg(28) = base_IMg * I_change;
    base_ab0 = pars_Mg(25);
    pars_Mg(25) = 2.7 * base_ab0;
        
    % ode15s
    [t,x] = ode15s(@(t,x) all_eqns_bp_Mg(t,x,pars_BP,pars_Mg,tchange,... 
                                                'ACEi',0, 'ARB', 0),...
                                                tspan,x1, options);                                      
    y=x(end,:);
    y_Mg = y(108:116);
    y_vals = y';
        
    % get Ca-Mg model fluxes
    flux_new = get_CaMg_fluxes(y_vals, y_Mg, pars_Mg);
        
    % saving results
    save_data_name = strcat('Rat_Data/rat_male_data_scenario_IMgi.mat');
    save(save_data_name, 'y_vals', 'flux_new')
end

%% 50% dietary calcium restriction simulation
if do_sim_Ca
    % read baseline datafile
    fname = './Rat_Data/rat_male_ss_data_scenario_normal_combined.mat';
    x1 = load(fname).SSdata;

    % parameters
    species = 'rat';
    sex = 'male';
    [pars_BP, pars_Mg, M] =  get_params_and_mass_matrix(x1,species,sex,htn_rsna,htn_raa,htn_renin,htn_ald);
        
    % options for ode
    options = odeset('Mass',M,'RelTol',1e-2,'AbsTol',1e-3*ones(1,length(x1)));

    % inhibiting dietary Ca intake by 70%
    base_ICa = pars_Mg(24);
    pars_Mg(24) = base_ICa * I_change;

    % ode15s
    [t,x] = ode15s(@(t,x) all_eqns_bp_Mg(t,x,pars_BP,pars_Mg,tchange,... 
                                                'ACEi',0, 'ARB', 0),...
                                                tspan,x1, options);                                      
    y=x(end,:);
    y_Mg = y(108:116);
    y_vals = y';
        
    % get Ca-Mg model fluxes
    flux_new = get_CaMg_fluxes(y_vals, y_Mg, pars_Mg);
        
    % saving results
    save_data_name = strcat('Rat_Data/rat_male_data_scenario_ICai.mat');
    save(save_data_name, 'y_vals', 'flux_new')
end

%% 25(OH)D inhibition simulation
if do_sim_vitD3
    % read baseline datafile
    fname = './Rat_Data/rat_male_ss_data_scenario_normal_combined.mat';
    x1 = load(fname).SSdata;    
        
    % parameters
    species = 'rat';
    sex = 'male';
    [pars_BP, pars_Mg, M] =  get_params_and_mass_matrix(x1,species,sex,htn_rsna,htn_raa,htn_renin,htn_ald);
        
    % options for ode
    options = odeset('Mass',M,'RelTol',1e-2,'AbsTol',1e-3*ones(1,length(x1)));

    % inhibiting dietary Mg intake by 70%
    base_vitD3 = pars_Mg(11);
    pars_Mg(11) = base_vitD3 * I_change;

    % ode15s
    [t,x] = ode15s(@(t,x) all_eqns_bp_Mg(t,x,pars_BP,pars_Mg,tchange,... 
                                                'ACEi',0, 'ARB', 0),...
                                                tspan,x1, options);                                      
    y=x(end,:);
    y_Mg = y(108:116);
    y_vals = y';
        
    % get Ca-Mg model fluxes
    flux_new = get_CaMg_fluxes(y_vals, y_Mg, pars_Mg);
        
    % saving results
    save_data_name = strcat('Rat_Data/rat_male_data_scenario_D3i.mat');
    save(save_data_name, 'y_vals', 'flux_new')
end


%% plot simulation results
if do_fig
    bpMg_indx = [51 37 50 64 80 83 1 6 9 33 109:112]; % BP=10, Mg=4; [50 63 79 82 1 5 7 8 32 59 108:111]
    
    % load the data
    fname_bl = './Rat_Data/rat_male_ss_data_scenario_normal_combined.mat';
    fname_Mg = strcat('Rat_Data/rat_male_data_scenario_IMgi.mat');
    fname_Ca = strcat('Rat_Data/rat_male_data_scenario_ICai.mat');
    fname_D3 = strcat('Rat_Data/rat_male_data_scenario_D3i.mat');
        
    x1_bl = load(fname_bl).SSdata;  
    x1_Mg = load(fname_Mg).y_vals;        
    x1_Ca = load(fname_Ca).y_vals;        
    x1_D3 = load(fname_D3).y_vals;
        
    x1_bl(109:112) = x1_bl(109:112)/(x1_bl(34)*1e-3);
    x1_Mg(109:112) = x1_Mg(109:112)/(x1_Mg(34)*1e-3);
    x1_Ca(109:112) = x1_Ca(109:112)/(x1_Ca(34)*1e-3);
    x1_D3(109:112) = x1_D3(109:112)/(x1_D3(34)*1e-3);
        
    x_pc_Mg = (x1_Mg(bpMg_indx) - x1_bl(bpMg_indx)) ./ x1_bl(bpMg_indx); 
    x_pc_Ca = (x1_Ca(bpMg_indx) - x1_bl(bpMg_indx)) ./ x1_bl(bpMg_indx);
    x_pc_D3 = (x1_D3(bpMg_indx) - x1_bl(bpMg_indx)) ./ x1_bl(bpMg_indx);

    %fractional changes in model variables
    x_pc_10 = zeros(14, 3);
    for ix=1:14
        x_pc_10(ix,:) = [x_pc_Mg(ix), x_pc_Ca(ix), x_pc_D3(ix)];
    end
    
    % Plot
    labels = {'','Mg^{2+} deficiency', 'Ca^{2+} deficiency', '25(OH)D deficiency'};
    
    x_label  = categorical(["MAP"; "CO"; "TPR"; "[ALD]"; "PRA"; "[AT1R-Ang II]"; "RSNA"; "RVR"; ...
                                "GFR"; "V_{ecf}";... 
                               "[PTH]"; "[1,25(OH)_2D_3]"; "[Mg^{2+}]"; "[Ca^{2+}]"]);

    mcolors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880;... 
        0.4940 0.1840 0.5560; 0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330; 0.6350 0.0780 0.1840; 0 0 1];
    
    graymap = gray(6);
    darkgray = graymap(2,:);
    
    lw = 3.0;
    f_gca = 16;
    fleg = 18;
    
    UP  = char(8593);
    DOWN  = char(8595);
    w = 0.9;

    t = tiledlayout(1,1,'TileSpacing','tight','Padding','Compact');
    ax1 = nexttile;
    hold on
    yline(0.0, 'color', darkgray, 'linewidth', 2.5)
    b = bar(x_pc_10, w, 'FaceColor','flat');
    b(1).CData = [0.5, 0, 0];
    b(2).CData = [0, 0, 1];
    b(3).CData = [0.0196, 0.4, 0.0314];
    set(gca, 'XTickLabel',x_label, 'XTick',1:numel(x_label), 'XTickLabelRotation', 0, 'fontsize', f_gca)
    legend(labels, 'fontsize', fleg, 'location','northeast')
    title('(A)', 'fontsize', 20)
    grid on
end

