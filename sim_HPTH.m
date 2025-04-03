%% list of tasks
sim_hpth = 0; % to simulate primary hyperparathyroididm
do_fig   = 1; % plot results

%% simulating hyperparathyroidism
if sim_hpth
    pth_incr = [2 3 5 7 10];

    % read baseline datafile
    IG = './Rat_Data/rat_male_ss_data_scenario_normal_combined.mat';
    x0 = load(IG).SSdata;
    flux_init = load(IG).fluxSS;
    
    species = 'rat';
    sex = 'male';
    [htn_rsna, htn_renin, htn_raa, htn_ald] = deal(1.3,1.0,1.0,1.0);
    [pars_BP, pars_Mg, M] =  get_params_and_mass_matrix(x0,species,sex,htn_rsna,htn_raa,htn_renin,htn_ald);
    options = odeset('Mass',M,'RelTol',1e-2,'AbsTol',1e-3*ones(1,length(x0)));

    for ip = 1:length(pth_incr)
        base_kprodPTH = pars_Mg(5);
        kprodPTH_new_val = pth_incr(ip) * base_kprodPTH;
        pars_Mg(5) = kprodPTH_new_val;
        
        tchange=0;
        tspan = [0 262800];
        [t,x] = ode15s(@(t,x) all_eqns_bp_Mg(t,x,pars_BP,pars_Mg,tchange,... 
                                                'ACEi',0, 'ARB', 0),...
                                                tspan,x0, options);
                                              
        y=x(end,:);
        y_Mg = y(108:116);
        y_vals = y';
        
        % get Ca-Mg model fluxes
        flux_new = get_CaMg_fluxes(y_vals, y_Mg, pars_Mg);

        pars_Mg(5) = base_kprodPTH;
    
        save_data_name = strcat('Rat_Data/rat_male_data_scenario_high_PTH_', num2str(pth_incr(ip)), '.mat');
        save(save_data_name, 'y_vals', 'flux_new')
    end  
end

if do_fig
    pth_incr = ["2" "3" "5" "7" "10"];
    
    flux_HP = containers.Map();
    x_HP_pc = containers.Map();

    x0(109:112) = x0(109:112)/(x0(34)*1e-3);
    
    bpMg_indx = [109:112 51 37 50 64 80 83 1 9 33]; % BP=9, Mg=4
    for ip = pth_incr
        fname = strcat('Rat_Data/rat_male_data_scenario_high_PTH_', ip, '.mat');
        x1 = load(fname).y_vals;
        x1(109:112) = x1(109:112)/(x1(34)*1e-3);
        
        flux1 = load(fname).flux_new;
        flux_HP(ip) = flux1;
        x1_pc = zeros(size(bpMg_indx));
        
        x_HP_pc(ip) = (x1(bpMg_indx) - x0(bpMg_indx)) ./ x0(bpMg_indx); % fractional change        
    end
    x_pc_10 = zeros(13, length(pth_incr));
    for ix=1:13
        x_temp = zeros(1,length(pth_incr));
        c = 0;
        for ip = pth_incr
            c=c+1;
            ax = x_HP_pc(ip);
            x_temp(c) = ax(ix);
        end
        x_pc_10(ix,:) = x_temp;
    end
    
    % Mg-Ca fluxes
    flux_pc_Mg = zeros(length(pth_incr)+1, 10);
    
    flux_pc_Mg(1,:) = [flux_init.Gut_absorption_Ca; flux_init.Gut_absorption_Mg;
                       flux_init.FastPool_to_Plasma_Ca; flux_init.FastPool_to_Plasma_Mg;
                       flux_init.Bone_resorption_Ca; flux_init.Bone_resorption_Mg;
                       flux_init.Plasma_to_FastPool_Ca; flux_init.Plasma_to_FastPool_Mg;
                       flux_init.Urine_excretion_Ca; flux_init.Urine_excretion_Mg]';
    
    for ii=2:length(pth_incr)+1
        flux_pc_Mg(ii,:) = [flux_HP(pth_incr(ii-1)).Gut_absorption_Ca;                                             %1
                            flux_HP(pth_incr(ii-1)).Gut_absorption_Mg;                                             %2
                            flux_HP(pth_incr(ii-1)).FastPool_to_Plasma_Ca;                                         %3
                            flux_HP(pth_incr(ii-1)).FastPool_to_Plasma_Mg;                                         %4
                            flux_HP(pth_incr(ii-1)).Bone_resorption_Ca;                                            %5
                            flux_HP(pth_incr(ii-1)).Bone_resorption_Mg;                                            %6
                            flux_HP(pth_incr(ii-1)).Plasma_to_FastPool_Ca;                                         %7
                            flux_HP(pth_incr(ii-1)).Plasma_to_FastPool_Mg;                                         %8
                            flux_HP(pth_incr(ii-1)).Urine_excretion_Ca;                                            %9
                            flux_HP(pth_incr(ii-1)).Urine_excretion_Mg]';                                          %10  
    end
    
    flux_pc_norm = flux_pc_Mg ./ flux_pc_Mg(1,:);
    
    % make figures
    labels = {'','PTH-2', 'PTH-3', 'PTH-5', 'PTH-7', 'PTH-10'};
    
    x_label  = categorical(["[PTH]"; "[1,25(OH)_2D_3]"; "[Mg^{2+}]"; "[Ca^{2+}]"; ...
                                "MAP"; "CO"; "TPR"; "[ALD]"; "PRA"; "[AT1R-Ang II]"; "RSNA"; ...
                                "GFR"; "V_{ecf}"]);
    
    % make figures
    mcolors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880;... 
        0.4940 0.1840 0.5560; 0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330; 0.6350 0.0780 0.1840; 0 0 1];
    
    graymap = gray(6);
    darkgray = graymap(2,:);
    
    lw = 3.0;
    %ls1 = '-';
    f_gca = 18;
    fleg = 18;
    
    UP  = char(8593);
    DOWN  = char(8595);
    w = 0.9;
    
    t = tiledlayout(2,3,'TileSpacing','Compact','Padding','Compact');
    ax1 = nexttile;
    hold on
    yline(0.0, 'color', darkgray, 'linewidth', 2.5)
    bar(x_pc_10(1:4,:), w)
    set(gca, 'XTickLabel',x_label(1:4), 'XTick',1:numel(x_label(1:4)), 'XTickLabelRotation', 0, 'fontsize', f_gca)
    ylim([0.0 5])
    legend(labels, 'fontsize', fleg, 'location','northeast')
    title('(A)', 'fontsize', 20)
    grid on

    ax2 = nexttile([1,2]);
    hold on
    yline(0.0, 'color', darkgray, 'linewidth', 2.5)
    bar(x_pc_10(5:end,:), w)
    set(gca, 'XTickLabel',x_label(5:end), 'XTick',1:numel(x_label(5:end)), 'XTickLabelRotation', 0, 'fontsize', f_gca)
    title('(B)', 'fontsize', 20)
    grid on
    
    ax3 = nexttile([1,3]);
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
    
    title('(C)', 'fontsize', 20)
    grid on
    
    ylabel(t,'Fractional change from baseline', 'fontsize', 20)
end