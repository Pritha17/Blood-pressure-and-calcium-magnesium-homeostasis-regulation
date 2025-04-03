%% hypertension scenarios
htn_sim         = 0;
%htn_type        = 'All';   % HTN-Combined
%htn_type        = 'RSNA';  % HTN-RSNA
%htn_type        = 'Renin'; % HTN-Renin
%htn_type        = 'ALD';   % HTN-ALD 
%htn_type        = 'RAA';   % HTN-AA

do_fig = 1;

% multiplicative factors for different types of HTN
if htn_sim
    [htn_rsna, htn_renin, htn_raa, htn_ald] = get_htn_factors(htn_type);
else
    [htn_rsna, htn_renin, htn_raa, htn_ald] = deal(1,1,1,1);
end

%% simulating hypertension
if htn_sim
    % read baseline datafile
    IG = './Rat_Data/rat_male_ss_data_scenario_normal_combined.mat';
    x0 = load(IG).SSdata;
    flux_init = load(IG).fluxSS;

    species = 'rat';
    sex = 'male';
    [pars_BP, pars_Mg, M] =  get_params_and_mass_matrix(x0,species,sex,htn_rsna,htn_raa,htn_renin,htn_ald);
    
    tchange=0;
    tspan = [0 600000];
    options = odeset('Mass',M,'RelTol',1e-2,'AbsTol',1e-3*ones(1,length(x0)));

    [t,x] = ode15s(@(t,x) all_eqns_bp_Mg(t,x,pars_BP,pars_Mg,tchange,... 
                                            'ACEi',0, 'ARB', 0),...
                                            tspan,x0, options);
                                          
    y=x(end,:);
    y_Mg = y(108:116);
    y_vals = y';   
    
    % get Ca-Mg model fluxes
    flux_new = get_CaMg_fluxes(y_vals, y_Mg, pars_Mg);

    save_data_name = strcat('Rat_Data/rat_male_data_scenario_HTN_', htn_type, '.mat');
    save(save_data_name, 'y_vals', 'flux_new')
end

if do_fig
    % read baseline datafile
    IG = './Rat_Data/rat_male_ss_data_scenario_normal_combined.mat';
    x0 = load(IG).SSdata;
    flux_init = load(IG).fluxSS;
    x0(109:112) = x0(109:112)/(x0(34)*1e-3);
    
    % read HTN datafiles
    htn_type = ["RSNA", "Renin", "ALD", "RAA", "All"];
    
    flux_htn = containers.Map();
    x_htn_pc = containers.Map();
    bpMg_indx = [51 64 80 83 1 6 8 9 33 60 109:112]; % BP=10, Mg=4

    for ht = htn_type
        fname = strcat('Rat_Data/rat_male_data_scenario_HTN_', ht, '.mat');
        x1 = load(fname).y_vals;
        x1(109:112) = x1(109:112)/(x1(34)*1e-3);
        flux_htn(ht) = load(fname).flux_new;
        x1_pc = zeros(size(bpMg_indx));
        x1_pc(:) = (x1(bpMg_indx) - x0(bpMg_indx)) ./ x0(bpMg_indx); % fractional change
        x_htn_pc(ht) = x1_pc;
    end
    
    % Model variables
    x_pc_10 = zeros(14, length(htn_type));
    for ix=1:14
        x_temp = zeros(1,length(htn_type));
        c = 0;
        for ht = htn_type
            c=c+1;
            ax = x_htn_pc(ht);
            x_temp(c) = ax(ix);
        end
        x_pc_10(ix,:) = x_temp;
    end
    
    % Mg-Ca fluxes
    flux_pc_Mg = zeros(length(htn_type)+1, 10);
    
    flux_pc_Mg(1,:) = [flux_init.Gut_absorption_Ca; flux_init.Gut_absorption_Mg;
                       flux_init.FastPool_to_Plasma_Ca; flux_init.FastPool_to_Plasma_Mg;
                       flux_init.Bone_resorption_Ca; flux_init.Bone_resorption_Mg;
                       flux_init.Plasma_to_FastPool_Ca; flux_init.Plasma_to_FastPool_Mg;
                       flux_init.Urine_excretion_Ca; flux_init.Urine_excretion_Mg]';
    
    for ii=2:length(htn_type)+1
        flux_pc_Mg(ii,:) = [flux_htn(htn_type(ii-1)).Gut_absorption_Ca;                                             %1
                            flux_htn(htn_type(ii-1)).Gut_absorption_Mg;                                             %2
                            flux_htn(htn_type(ii-1)).FastPool_to_Plasma_Ca;                                         %3
                            flux_htn(htn_type(ii-1)).FastPool_to_Plasma_Mg;                                         %4
                            flux_htn(htn_type(ii-1)).Bone_resorption_Ca;                                            %5
                            flux_htn(htn_type(ii-1)).Bone_resorption_Mg;                                            %6
                            flux_htn(htn_type(ii-1)).Plasma_to_FastPool_Ca;                                         %7
                            flux_htn(htn_type(ii-1)).Plasma_to_FastPool_Mg;                                         %8
                            flux_htn(htn_type(ii-1)).Urine_excretion_Ca;                                            %9
                            flux_htn(htn_type(ii-1)).Urine_excretion_Mg]';                                          %10  
    end
    
    flux_pc_norm = flux_pc_Mg ./ flux_pc_Mg(1,:);
    
    % Plot
    labels = {'', 'HTN-RSNA', 'HTN-Renin', 'HTN-ALD', 'HTN-AA', 'HTN-Combined'};
    
    x_label  = categorical(["MAP"; "[ALD]"; "PRA"; "[AT1R-Ang II]"; "RSNA"; "RVR"; ...
                                "RBF"; "GFR"; "V_{ecf}"; "T_{Na}"; ... 
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

    t = tiledlayout(3,1,'TileSpacing','tight','Padding','Compact');
    ax1 = nexttile;
    hold on
    yline(0.0, 'color', darkgray, 'linewidth', 2.5)
    bar(x_pc_10, w, 'FaceColor','flat'); 
    set(gca,'TickLength',[0 0], 'fontsize', f_gca) % removes tick marks and 
    set(gca,'Xtick',[])
    set(gca, 'YGrid', 'on', 'XGrid', 'off')
    ylim([5 11])
    legend(labels, 'fontsize', fleg, 'location','northeast')
    title('(A)', 'fontsize', 20)

    ax2 = nexttile(2,[(3-2) 1]);
    hold on
    yline(0.0, 'color', darkgray, 'linewidth', 2.5)
    bar(x_pc_10, w, 'FaceColor','flat');
    set(gca, 'XTickLabel',x_label, 'XTick',1:numel(x_label), 'XTickLabelRotation', 0, 'fontsize', f_gca)
    set(gca, 'YGrid', 'on', 'XGrid', 'off')
    ylim([-0.5 1.0])
    
    ax3 = nexttile;
    hold on
    yline(0.0, 'color', darkgray, 'linewidth', 2.5)
    x_SS_2 = categorical(["Gut absorption (Ca^{2+})"; "Gut absorption (Mg^{2+})";...
                        "Bone to plasma (Ca^{2+})"; "Bone to plasma (Mg^{2+})";...
                        "Bone resorption (Ca^{2+})"; "Bone resorption (Mg^{2+})";...
                        "Plasma to bone (Ca^{2+})"; "Plasma to bone (Mg^{2+})";...
                         "Urine excretion (Ca^{2+})"; "Urine excretion (Mg^{2+})"]);
    
    y_SS_2 = [flux_pc_norm(2:end,1).'-1; flux_pc_norm(2:end,2).'-1; 
              flux_pc_norm(2:end,3).'-1; flux_pc_norm(2:end,4).'-1; 
              flux_pc_norm(2:end,5).'-1; flux_pc_norm(2:end,6).'-1;
              flux_pc_norm(2:end,7).'-1; flux_pc_norm(2:end,8).'-1; 
              flux_pc_norm(2:end,9).'-1; flux_pc_norm(2:end,10).'-1];
    bar(y_SS_2, w)
    set(gca, 'XTickLabel',x_SS_2, 'XTick',1:numel(x_SS_2), 'XTickLabelRotation', 45, 'fontsize', f_gca)
    set(gca, 'YGrid', 'on', 'XGrid', 'off')
    title('(B)', 'fontsize', 20)
    
    ylabel(t,'Fractional change from baseline', 'fontsize', 20)
end

