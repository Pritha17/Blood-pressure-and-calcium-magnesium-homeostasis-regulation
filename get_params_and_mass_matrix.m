function [pars_BP, pars_Mg, M] =  get_params_and_mass_matrix(x1,species,sex,htn_rsna,htn_raa,htn_renin,htn_ald)
        pars_BP = get_pars(species,sex,'RSNA',htn_rsna,'R_aa',htn_raa,'Renin',htn_renin,'ALD',htn_ald);
        run('read_params.m')
        params_Mg = param_vals;
        pars_Mg = params_Mg(:); % converting to column vector
    
        % ensure consistency
        total_PTH   = x1(109);
        total_D3    = x1(110);
        total_Mg    = x1(111);
        total_Ca    = x1(112);
        V_b         = x1(34)*1e-3;
        plasmaPTH   = total_PTH/V_b;
        plasmaD3    = total_D3/V_b;
        plasmaMg    = total_Mg/V_b;
        plasmaCa    = total_Ca/V_b;
        pars_BP(44) = plasmaCa;
        pars_BP(45) = plasmaMg;
        pars_BP(46) = plasmaD3;
        pars_BP(47) = plasmaPTH;
        
        pars_Mg(34) = x1(9)*1e-3;  % GFR_base
        pars_Mg(4)  = x1(34)*1e-3; % plasma volume
        pars_Mg(76) = x1(64);      % C_ALD
        pars_Mg(77) = x1(17);      % gamma_at
        pars_Mg(78) = x1(24);      % psi_al
        pars_Mg(80) = x1(83);      % AT1R-bound AngII
        
        % ODE options
        % defining the mass matrix
        M = zeros(length(x1));
        ind = [33 39 55 57 59 60 65 74 79 81:86 108:116];
        for i=1:length(ind)
            M(ind(i),ind(i)) = 1; 
        end
        % now fill in extra entries
        M(55,53) = -3/4;  % a_baro, a_auto
        M(59,38) = -0.2;  % delta_ra, P_ra
end