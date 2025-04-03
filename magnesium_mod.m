function dydt = magnesium_mod(t,y,params)

% get variable inputs
Ca_fixed          = false; % turn on to have no calcium dynamics
Mg_fixed          = false; % turn on to have no magnesium dynamics
PTHg_fixed        = false; % turn on to have no PTHg dynamics
PTHp_fixed        = false; % turn on to have no PTHp dynamics
D3_fixed          = false; % turn on to have no D3 dynamics
NCaf_fixed        = false; % turn on to have no NCaf dynamics
NCas_fixed        = true; % turn on to have no NCas dynamics
NMgf_fixed        = false; % turn on to have no NMgf dynamics
NMgs_fixed        = true; % turn on to have no NMgs dynamics

% variable names
PTH_g         = y(1); % amount of PTH in parathyroid gland pool (pmol)
PTH_p         = y(2); % amount of PTH in plasma (pmol)
D3_p          = y(3); % amount of vitamin D3 in plasma (pmol)
Mg_p          = y(4); % amount of magnesium in plasma (mmol)
Ca_p          = y(5); % amount of calcium in plasma (mmol)
NCa_f         = y(6); % amount of calcium in the rapidly excangeable pool (mmol)
NCa_s         = y(7); % amount of calcium in the slow pool (mmol)
NMg_f         = y(8); % amount of magnesium in the rapidly excangeable pool (mmol)
NMg_s         = y(9); % amount of magnesium in the slow pool (mmol)

% parameters
gamma_p = params(1);
beta_exo_PTHg = params(2);
Cm_half = params(3);
Vp = params(4);
k_prod_PTHg = params(5);
gamma_Ca = params(6);
gamma_Mg = params(7);
gamma_prod_D3 = params(8);
k_PTHg_deg = params(9);
k_PTHp_deg = params(10);
D3_inact = params(11);
K_D3 = params(12);
gamma_deg_PTHp = params(13);
k_deg_D3 = params(14);
nconv = params(15);
K_conv_PTH = params(16);
gamma_conv_Ca = params(17);
gamma_conv_D3 = params(18);
delta_Mg_act = params(19);
K1_Mg_act = params(20);
K2_Mg_act = params(21);
k_min_conv = params(22);
delta_max_D3act = params(23);
ICa = params(24);
Gamma_abs0 = params(25);
delta_abs_D3 = params(26);
K_abs_D3 = params(27);
IMg = params(28);
gamma_passive_Mg = params(29);
gamma_active_Mg = params(30);
gamma_D3_Mg = params(31);
V_active = params(32);
K_active = params(33);
GFR = params(34);
Lambda_PT0_Ca = params(35);
delta_PT_max_Ca = params(36); 
Lambda_TAL0_Ca = params(37); 
delta_TAL_max_Ca = params(38);  
Lambda_DCT0_Ca = params(39);
delta_DCT_max_Ca = params(40); 
Lambda_PT0_Mg = params(41);
delta_PT_max_Mg = params(42); 
Lambda_TAL0_Mg = params(43); 
delta_TAL_max_Mg = params(44);
Lambda_DCT0_Mg = params(45);
delta_DCT_max_Mg = params(46);  
PTHp_ref = params(47); 
K_TAL_PTHp = params(48); 
K_DCT_PTHp = params(49);
K_DCT_D3p  = params(50);
nPT = params(51);
nTAL = params(52);
Cap_ref = params(53);
Mgp_ref = params(54);
Gamma_res_min = params(55);
delta_res_max = params(56);
K_PTHp_res = params(57);
K_D3p_res = params(58);
kappa_b_Ca = params(59);
kappa_b_Mg = params(60);
k_pf_Ca = params(61);
k_fp_Ca = params(62);
Gamma_ac_Ca = params(63);
k_pf_Mg = params(64);
k_fp_Mg = params(65);
Gamma_ac_Mg = params(66);
PTHg_norm = params(67);
PTHp_norm = params(68);
D3p_norm = params(69);
Mgp_norm = params(70);
Cap_norm = params(71);
NCaf_norm = params(72);
NCas_norm = params(73);
NMgf_norm = params(74);
NMgs_norm = params(75);
C_al = params(76);
gamma_at = params(77);
psi_al = params(78);
AT1R_eq = params(79);
AT1R = params(80);
K_PTH_ALD = params(81);
K_PTH_AT1R = params(82);

dydt = zeros(length(y),1);

% change to concentrations
PTHp_con = PTH_p/Vp;
Cap_con = Ca_p/Vp;
Mgp_con = Mg_p/Vp;
D3p_con  = D3_p/Vp;

% All equations have been normalized
% Parathyroid glands
if PTHg_fixed
    dydt(1) = 0;
else
    h = (Mgp_con/0.4)^50 / (1 + (Mgp_con/0.4)^50);

    PTHg_prod_basal = k_prod_PTHg / PTHg_norm;
    PTHg_prod_effect_D3 = 1 / (1 + gamma_prod_D3 * D3p_norm * D3p_con);
    PTHg_synthesis = PTHg_prod_basal * PTHg_prod_effect_D3;

    PTHg_deg = k_PTHg_deg * PTH_g;

    normal_high_Mg_effect = ((gamma_Ca / (Cap_norm*Cap_con)) ^ gamma_p) * ...
    (beta_exo_PTHg - 1 / (1 + (Cm_half/Mgp_con)^2));
    
    low_Mg_effect = gamma_Mg / (1 + 0.25/Mgp_con);
    
    F_Ca_Mg = h * normal_high_Mg_effect + (1-h) * low_Mg_effect;
    
    MRB = 0;
    ARB = 0;
    F_ALD = (1-MRB) * (1.72 * C_al^8 / (K_PTH_ALD^8 + C_al^8)); % MRB: mineralocorticoid blocker
    F_AT1R = (1-ARB) * (1.67 * AT1R^8 / (K_PTH_AT1R^8 + AT1R^8));
    F_ALD_ATIR = min(2.7, 1 + F_ALD + F_AT1R);

    PTHg_exocytosis = F_ALD_ATIR * F_Ca_Mg * PTH_g;

    dydt(1) = PTHg_synthesis - PTHg_exocytosis - PTHg_deg;
end

% Plasma PTH
if PTHp_fixed
    dydt(2) = 0;
else
    PTHp_deg = k_PTHp_deg * PTHp_con;
    PTHp_influx = PTHg_exocytosis * PTHg_norm / PTHp_norm;

    dydt(2) = PTHp_influx - PTHp_deg;
end

% Vitamin D3
if D3_fixed
    dydt(3) = 0;
else
    PTH_impact_D3_act = PTHp_con^nconv / (PTHp_con^nconv + (K_conv_PTH/PTHp_norm)^nconv);
    Ca_impact_D3_act  = 1 / (1 + gamma_conv_Ca * Cap_norm * Cap_con);
    D3_impact_D3_act  = 1 / (1 + gamma_conv_D3 * D3p_norm * D3p_con);
    Mg_impact_D3_act_1 = delta_Mg_act * Mgp_con^4 / ((K1_Mg_act/Mgp_norm)^4 + Mgp_con^4);
    Mg_impact_D3_act_2 = delta_Mg_act * (K2_Mg_act/Mgp_norm)^4 / ((K2_Mg_act/Mgp_norm)^4 + Mgp_con^4);
    h_m = 1 / (1 + (Mgp_con/2.5)^50);
    
    Mg_impact_D3_act = h_m * Mg_impact_D3_act_1 + (1-h_m) * Mg_impact_D3_act_2;
    
    F_D3_act = PTH_impact_D3_act * Ca_impact_D3_act * D3_impact_D3_act * Mg_impact_D3_act; %
    
    D3_influx = (k_min_conv + delta_max_D3act * F_D3_act) * (D3_inact/D3p_norm);
    
    PTH_impact_D3_deg = 1 / (1 + gamma_deg_PTHp * PTHp_norm * PTHp_con);
    Mg_impact_D3_inact = Mgp_con^4 / ((K_D3/Mgp_norm)^4 + Mgp_con^4);
     
    deg_D3_act = (k_deg_D3 * (PTH_impact_D3_deg + Mg_impact_D3_inact)) * D3p_con;
    
    dydt(3) = D3_influx - deg_D3_act;
end

% Intestinal compartment
% Effect of vitamin D3
Gut_impact_D3 = D3p_con^2/((K_abs_D3/D3p_norm)^2 + D3p_con^2);

% intestinal calcium absorption
Gut_frac_absorption_Ca = Gamma_abs0 + delta_abs_D3*Gut_impact_D3;
Gut_absorption_Ca = (ICa*Gut_frac_absorption_Ca)/Cap_norm;

% intestinal magnesium absorption
passive_Mg_absorption = gamma_passive_Mg;
active_Mg_absorption = gamma_active_Mg * V_active / (K_active + IMg);
D3_effect_Mg_absorption = gamma_D3_Mg * Gut_impact_D3;
Gut_frac_absorption_Mg = passive_Mg_absorption + active_Mg_absorption + D3_effect_Mg_absorption;
Gut_absorption_Mg = (IMg * Gut_frac_absorption_Mg) / Mgp_norm;

% Kidney compartment
% Effect of PTH in PT
delta_PT_PTH = 1/(1 + (PTHp_con * PTHp_norm/PTHp_ref)^nPT);

% Effect of CaSR and PTH in TAL
delta_TAL_CaSR = 1/((1 + (Cap_con * Cap_norm/Cap_ref)^nTAL)*(1 + 0.6*(Mgp_con * Mgp_norm/Mgp_ref)^nTAL));
delta_TAL_PTH = PTHp_con/(PTHp_con + K_TAL_PTHp/PTHp_norm);

% Effect of PTH and vitamin D3 in DCT
delta_DCT_PTH = PTHp_con/(PTHp_con + K_DCT_PTHp/PTHp_norm);
delta_DCT_D3  = D3p_con/(D3p_con + K_DCT_D3p/D3p_norm);

% renal calcium handling
Renal_filtration_Ca = GFR * Cap_con;
Lambda_PT_Ca = Lambda_PT0_Ca + delta_PT_max_Ca * (0.5*delta_PT_PTH + 0.5*gamma_at);
Lambda_TAL_Ca = Lambda_TAL0_Ca + delta_TAL_max_Ca * (0.55*delta_TAL_CaSR + 0.25*delta_TAL_PTH + 0.2*gamma_at);
Lambda_DCT_Ca = Lambda_DCT0_Ca + delta_DCT_max_Ca * (0.6*delta_DCT_PTH + 0.2*delta_DCT_D3 + 0.2*psi_al);
Renal_frac_reab_Ca = min(0.995, Lambda_PT_Ca + Lambda_TAL_Ca + Lambda_DCT_Ca);
Urine_excretion_Ca = (1-Renal_frac_reab_Ca)*Renal_filtration_Ca;

% renal magnesium handling
Renal_filtration_Mg = GFR * Mgp_con;
Lambda_PT_Mg = Lambda_PT0_Mg + delta_PT_max_Mg * (0.5*delta_PT_PTH + 0.5*gamma_at);
Lambda_TAL_Mg = Lambda_TAL0_Mg + delta_TAL_max_Mg * (0.55*delta_TAL_CaSR + 0.25*delta_TAL_PTH + 0.2*gamma_at);
Lambda_DCT_Mg = Lambda_DCT0_Mg + delta_DCT_max_Mg*(0.6*delta_DCT_PTH + 0.2*delta_DCT_D3 + 0.2*psi_al);
Renal_frac_reab_Mg = min(0.995, Lambda_PT_Mg + Lambda_TAL_Mg + Lambda_DCT_Mg);
Urine_excretion_Mg = (1-Renal_frac_reab_Mg)*Renal_filtration_Mg;

% Bone compartment
% bone calcium and magnesium resorption
PTHp_res_effect = delta_res_max*0.2*PTHp_con^2/(PTHp_con^2 + (K_PTHp_res/PTHp_norm)^2);
D3_res_effect   = delta_res_max*0.6*D3p_con^2/(D3p_con^2 + (K_D3p_res/D3p_norm)^2);
AngII_res_effect = delta_res_max*0.2*AT1R^2/(AT1R^2 + (AT1R_eq)^2);
Bone_resorption_Ca = Gamma_res_min + PTHp_res_effect + D3_res_effect + AngII_res_effect;
Bone_resorption_Mg = 0.025 * (Gamma_res_min + PTHp_res_effect + D3_res_effect + AngII_res_effect);

% bone fast pool calcium
Plasma_to_FastPool_Ca = k_pf_Ca*Cap_con*Vp;
FastPool_to_Plasma_Ca = k_fp_Ca*NCa_f;
Bone_accretion_Ca = Gamma_ac_Ca*NCa_f;

% bone fast pool magnesium
Plasma_to_FastPool_Mg = k_pf_Mg*Mgp_con*Vp;
FastPool_to_Plasma_Mg = k_fp_Mg*NMg_f;
Bone_accretion_Mg = Gamma_ac_Mg*NMg_f;

% plasma magnesium
if Mg_fixed
    dydt(4) = 0;
else
    dydt(4) = (1 - kappa_b_Mg)*(Gut_absorption_Mg + Bone_resorption_Mg/Mgp_norm + ...
        FastPool_to_Plasma_Mg - Plasma_to_FastPool_Mg - Urine_excretion_Mg) / Vp;
end

% plasma calcium
if Ca_fixed
    dydt(5) = 0;
else
    dydt(5) = (1 - kappa_b_Ca)*(Gut_absorption_Ca + Bone_resorption_Ca/Cap_norm + ...
        FastPool_to_Plasma_Ca - Plasma_to_FastPool_Ca - Urine_excretion_Ca) / Vp;
end

% rapidly exchangeable calcium pool
if NCaf_fixed
    dydt(6) = 0;
else
    dydt(6) = Plasma_to_FastPool_Ca - FastPool_to_Plasma_Ca - Bone_accretion_Ca;
end

% slowly exchangeable calcium pool
if NCas_fixed
    dydt(7) = 0;
else
    dydt(7) = Bone_accretion_Ca - Bone_resorption_Ca/NCas_norm;
end

% rapidly exchangeable magnesium pool
if NMgf_fixed
    dydt(8) = 0;
else
    dydt(8) = Plasma_to_FastPool_Mg - FastPool_to_Plasma_Mg - Bone_accretion_Mg;
end

% slowly exchangeable magnesium pool
if NMgs_fixed
    dydt(9) = 0;
else
    dydt(9) = Bone_accretion_Mg - Bone_resorption_Mg/NMgs_norm;
end

end