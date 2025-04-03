function flux_new = get_CaMg_fluxes(y_vals, y_Mg, pars_Mg)
        % ensure consistency in parameters
        pars_Mg(34) = y_vals(9)*1e-3;  % GFR_base
        pars_Mg(4)  = y_vals(34)*1e-3; % plasma volume
        pars_Mg(76) = y_vals(64);      % C_ALD
        pars_Mg(77) = y_vals(17);      % gamma_at
        pars_Mg(78) = y_vals(24);      % psi_al
        pars_Mg(80) = y_vals(83);      % AT1R-bound AngII
        
        flux_new = compute_fluxes(y_Mg, pars_Mg);
end