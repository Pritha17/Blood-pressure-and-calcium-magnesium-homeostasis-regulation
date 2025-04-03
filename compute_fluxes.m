function vals = compute_fluxes(yvals,params)

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

% variable names
PTH_g         = yvals(:,1); % amount of PTH in parathyroid gland pool (pmol)
PTH_p         = yvals(:,2); % amount of PTH in plasma (pmol)
D3_p          = yvals(:,3); % amount of vitamin D3 in plasma (pmol)
Mg_p          = yvals(:,4); % amount of magnesium in plasma (mmol)
Ca_p          = yvals(:,5); % amount of calcium in plasma (mmol)
NCa_f         = yvals(:,6); % amount of calcium in the rapidly excangeable pool (mmol)
NCa_s         = yvals(:,7); % amount of calcium in the slow pool (mmol)
NMg_f         = yvals(:,8); % amount of magnesium in the rapidly excangeable pool (mmol)
NMg_s         = yvals(:,9); % amount of magnesium in the slow pool (mmol)

% change to concentrations
PTHp_con = PTH_p./Vp;
Cap_con = Ca_p./Vp;
Mgp_con = Mg_p./Vp;
D3p_con  = D3_p./Vp;

% Parathyroid glands
h = (Mgp_con/0.4)^50 / (1 + (Mgp_con/0.4)^50);

vals.PTHg_prod_basal = k_prod_PTHg / PTHg_norm;
vals.PTHg_prod_effect_D3 = 1 ./ (1 + gamma_prod_D3 * D3p_norm * D3p_con);
vals.PTHg_synthesis = vals.PTHg_prod_basal * vals.PTHg_prod_effect_D3;

vals.PTHg_deg = k_PTHg_deg * PTH_g;

vals.normal_high_Mg_effect = ((gamma_Ca ./ (Cap_norm*Cap_con)) ^ gamma_p) * ...
    (beta_exo_PTHg - 1 ./ (1 + (Cm_half/Mgp_con)^2));

vals.low_Mg_effect = gamma_Mg ./ (1 + 0.25/Mgp_con);
vals.F_Ca_Mg = h * vals.normal_high_Mg_effect + (1-h) * vals.low_Mg_effect;
MRB=0;
ARB = 0;
vals.F_ALD = (1-MRB) * (1.72 * C_al^8 / (K_PTH_ALD^8 + C_al^8));
vals.F_AT1R = (1-ARB) * (1.67 * AT1R^8 / (K_PTH_AT1R^8 + AT1R^8));
vals.F_ALD_ATIR = min(2.7, 1 + vals.F_ALD + vals.F_AT1R);
vals.PTHg_exocytosis = vals.F_ALD_ATIR * vals.F_Ca_Mg * PTH_g;

% Plasma PTH
vals.PTHp_deg = k_PTHp_deg * PTHp_con;
vals.PTHp_influx = vals.PTHg_exocytosis * PTHg_norm / PTHp_norm;

% Vitamin D3
vals.PTH_impact_D3_act = PTHp_con.^nconv ./ (PTHp_con.^nconv + (K_conv_PTH/PTHp_norm).^nconv);
vals.Ca_impact_D3_act  = 1 ./ (1 + gamma_conv_Ca * Cap_norm * Cap_con);
vals.D3_impact_D3_act  = 1 ./ (1 + gamma_conv_D3 * D3p_norm * D3p_con);
vals.Mg_impact_D3_act_1 = delta_Mg_act * Mgp_con.^4 / ((K1_Mg_act/Mgp_norm).^4 + Mgp_con.^4);
vals.Mg_impact_D3_act_2 = delta_Mg_act * (K2_Mg_act/Mgp_norm)^4 / ((K2_Mg_act/Mgp_norm)^4 + Mgp_con^4);
h_m = 1 / (1 + (Mgp_con/2.5)^50);
vals.Mg_impact_D3_act = h_m * vals.Mg_impact_D3_act_1 + (1-h_m) * vals.Mg_impact_D3_act_2;

vals.F_D3_act = vals.PTH_impact_D3_act * vals.Ca_impact_D3_act * vals.D3_impact_D3_act * vals.Mg_impact_D3_act; % unitless
vals.D3_influx = (k_min_conv + delta_max_D3act * vals.F_D3_act) * (D3_inact/D3p_norm);

vals.PTH_impact_D3_deg = 1 ./ (1 + gamma_deg_PTHp * PTHp_norm * PTHp_con);
vals.Mg_impact_D3_inact = Mgp_con.^4 / ((K_D3/Mgp_norm).^4 + Mgp_con.^4);
vals.deg_D3_act = (k_deg_D3 * (vals.PTH_impact_D3_deg + vals.Mg_impact_D3_inact)) * D3p_con;

% Intestinal compartment
% Effect of vitamin D3
vals.Gut_impact_D3 = D3p_con.^2 ./ ((K_abs_D3/D3p_norm).^2 + D3p_con.^2);

% intestinal calcium absorption
vals.Gut_frac_absorption_Ca = Gamma_abs0 + delta_abs_D3 * vals.Gut_impact_D3;
vals.Gut_absorption_Ca = (ICa * vals.Gut_frac_absorption_Ca) / Cap_norm;

% intestinal magnesium absorption
vals.passive_Mg_absorption = gamma_passive_Mg;
vals.active_Mg_absorption = gamma_active_Mg * V_active / (K_active + IMg);
vals.D3_effect_Mg_absorption = gamma_D3_Mg * vals.Gut_impact_D3;
vals.Gut_frac_absorption_Mg = vals.passive_Mg_absorption + vals.active_Mg_absorption + vals.D3_effect_Mg_absorption;
vals.Gut_absorption_Mg = (IMg * vals.Gut_frac_absorption_Mg) / Mgp_norm;

% Kidney compartment
% Effect of PTH in PT
vals.delta_PT_PTH = 1 ./ (1 + (PTHp_con * PTHp_norm / PTHp_ref).^nPT);

% Effect of CaSR and PTH in TAL
vals.delta_TAL_CaSR = 1 ./ ((1 + (Cap_con * Cap_norm / Cap_ref).^nTAL) * (1 + 0.6 * (Mgp_con * Mgp_norm / Mgp_ref).^nTAL));
vals.delta_TAL_PTH = PTHp_con ./ (PTHp_con + K_TAL_PTHp/PTHp_norm);

% Effect of PTH and vitamin D3 in DCT
vals.delta_DCT_PTH = PTHp_con ./ (PTHp_con + K_DCT_PTHp/PTHp_norm);
vals.delta_DCT_D3  = D3p_con ./ (D3p_con + K_DCT_D3p/D3p_norm);

% renal calcium handling
%disp(Cap_con)
vals.Renal_filtration_Ca = GFR * Cap_con;
vals.Lambda_PT_Ca = Lambda_PT0_Ca + delta_PT_max_Ca * (0.5*vals.delta_PT_PTH + 0.5*gamma_at);
vals.Lambda_TAL_Ca = Lambda_TAL0_Ca + delta_TAL_max_Ca * (0.55*vals.delta_TAL_CaSR + 0.25*vals.delta_TAL_PTH + 0.2*gamma_at);
vals.Lambda_DCT_Ca = Lambda_DCT0_Ca + delta_DCT_max_Ca * (0.6*vals.delta_DCT_PTH + 0.2*vals.delta_DCT_D3 + 0.2*psi_al);
vals.Renal_frac_reab_Ca = min(0.995, vals.Lambda_PT_Ca + vals.Lambda_TAL_Ca + vals.Lambda_DCT_Ca);
vals.Urine_excretion_Ca = (1 - vals.Renal_frac_reab_Ca) * vals.Renal_filtration_Ca;

% renal magnesium handling
vals.Renal_filtration_Mg = GFR * Mgp_con;
vals.Lambda_PT_Mg = Lambda_PT0_Mg + delta_PT_max_Mg * (0.5*vals.delta_PT_PTH + 0.5*gamma_at);
vals.Lambda_TAL_Mg = Lambda_TAL0_Mg + delta_TAL_max_Mg * (0.55*vals.delta_TAL_CaSR + 0.25*vals.delta_TAL_PTH + 0.2*gamma_at);
vals.Lambda_DCT_Mg = Lambda_DCT0_Mg + delta_DCT_max_Mg*(0.6*vals.delta_DCT_PTH + 0.2*vals.delta_DCT_D3 + 0.2*psi_al);
vals.Renal_frac_reab_Mg = min(0.995, vals.Lambda_PT_Mg + vals.Lambda_TAL_Mg + vals.Lambda_DCT_Mg);
vals.Urine_excretion_Mg = (1 - vals.Renal_frac_reab_Mg) * vals.Renal_filtration_Mg;

% Bone compartment
% bone calcium resorption
vals.PTHp_res_effect = delta_res_max*0.2*PTHp_con.^2./(PTHp_con.^2 + (K_PTHp_res/PTHp_norm).^2);
vals.D3_res_effect   = delta_res_max*0.6*D3p_con.^2./(D3p_con.^2 + (K_D3p_res/D3p_norm).^2);
vals.AngII_res_effect = delta_res_max*0.2*AT1R^2/(AT1R^2 + AT1R_eq^2);
vals.Bone_resorption_Ca = Gamma_res_min + vals.PTHp_res_effect + vals.D3_res_effect + vals.AngII_res_effect;
if NMg_s <= 0
    vals.Bone_resorption_Mg = 0;
else
    vals.Bone_resorption_Mg = 0.025 * (Gamma_res_min + vals.PTHp_res_effect + vals.D3_res_effect + vals.AngII_res_effect);
end

% bone fast pool calcium
vals.Plasma_to_FastPool_Ca = k_pf_Ca*Cap_con*Vp;
vals.FastPool_to_Plasma_Ca = k_fp_Ca*NCa_f;
vals.Bone_accretion_Ca = Gamma_ac_Ca*NCa_f;

% bone fast pool magnesium
vals.Plasma_to_FastPool_Mg = k_pf_Mg*Mgp_con*Vp;
vals.FastPool_to_Plasma_Mg = k_fp_Mg*NMg_f;
vals.Bone_accretion_Mg = Gamma_ac_Mg*NMg_f;

end