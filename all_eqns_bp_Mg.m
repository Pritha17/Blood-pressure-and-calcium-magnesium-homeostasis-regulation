function f = all_eqns_bp_Mg(t,x,pars_BP,pars_Mg,tchange,varargin)
    num_BP_var = 107;  % number of variables in BP model;
    x_BP = x(1:num_BP_var);
    x_Mg  = x(num_BP_var+1:end);

    % set plasma [PTH], [D3], [Ca], and [Mg] in pars_BP
    total_PTH   = x(num_BP_var+2);
    total_D3    = x(num_BP_var+3);
    total_Mg    = x(num_BP_var+4);
    total_Ca    = x(num_BP_var+5);
    V_b         = x(34)*1e-3;
    plasmaPTH   = total_PTH/V_b;
    plasmaD3    = total_D3/V_b;
    plasmaMg    = total_Mg/V_b;
    plasmaCa    = total_Ca/V_b;
    pars_BP(44) = plasmaCa;
    pars_BP(45) = plasmaMg;
    pars_BP(46) = plasmaD3;
    pars_BP(47) = plasmaPTH;
    x_BP0 = zeros(num_BP_var,1);
    
    f_BP = Ca_Mg_bp_reg_mod(t,x_BP,x_BP0,pars_BP,tchange,varargin{:});
    
    % set GFR and plasma volume in pars_Mg
    pars_Mg(34) = x(9)*1e-3;  % GFR_base
    pars_Mg(4)  = x(34)*1e-3; % plasma volume
    pars_Mg(76) = x(64);      % C_ALD
    pars_Mg(77) = x(17);      % gamma_at
    pars_Mg(78) = x(24);      % psi_al
    pars_Mg(80) = x(83);      % AT1R-bound AngII
    
    f_Mg  = magnesium_mod(t,x_Mg,pars_Mg);
    
    f = [f_BP; f_Mg];
end