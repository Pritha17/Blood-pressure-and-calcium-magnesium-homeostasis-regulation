% This is a model of long-term model blood pressure regulation.
% It is adopted with modifications from Karaaslan 2005, Leete 2018, and Ahmed 2020.

% This function file is to solve for the steady state solution in different scenarios.

% Input
% t        - time
% x        - variables
% x_p      - variable derivatives
% pars     - parameters
% tchange  - time at which to change RPP in simulation
% varargin - cell containing different scenarios and the corresponding parameters 
%            it is a variable length input argument

% Output
% f        - left hand side of f(t,x(t),x'(t);theta) = 0.

function f = Ca_Mg_bp_reg_mod(t,x,x_p,pars,tchange,varargin)

%% Retrieve species and sex identifier. 
spe_par = pars(1);
sex_par = pars(2);

if     spe_par == 1
    species = 'human';
elseif spe_par == 0
    species = 'rat';
end

if     sex_par == 1
    sex = 'male';
elseif sex_par == 0
    sex = 'female';
end

%% Retrieve species specific additional things.

if strcmp(species, 'rat')
    % Steady state variable values
    SSdata_input     = pars(end-106:end); % 106 = no. of variables-1
    pars(end-106:end) = '';
    
    % Fixed variable parameters
    fixed_var_pars  = pars(end-8:end);
    pars(end-8:end) = '';
end

% Scaling factors
% Rat value = Human value x SF
% Note: This includes conversion of units.
SF_S = pars(end-3); % sodium flow
SF_U = pars(end-2); % urine flow
SF_R = pars(end-1); % resistance
SF_V = pars(end  ); % volume

%% Default parameter inputs for changes in simulation.

% Drugs; variable names changed by Pritha
kappa_AngII = 0; % Ang II infusion rate fmol/ml/min
kappa_ACEi  = 0; % ACE inhibitor %
kappa_ARB   = 0; % Angiotensin receptor blocker %
kappa_f     = 0; % Furosemide values. array of length 2
kappa_f_md  = 0; % Furosemide values. array of length 2
NSAID           = 0; % NSAID indicator. O for none, 1 for normal, 2 for high dose.

% Renal perfusion pressure
RPP_ind = false; % Boolean for controlling renal perfusion pressure
RPP_per = 0    ; % Renal perfusion pressure mmHg

% Autoregulatory mechanisms
denerve     = false; % Boolean for unilateral renal denervation by fixing rsna = 1, which is baseline.
no_TGF      = false; % Boolean for canceling the tubuloglomerular feedback.
no_Myo      = false; % Boolean for canceling the myogenic response.
lin_Myo     = false; % Boolean for having linear myogenic response.
imp_Myo_ind = false; % Boolean for having impaired myogenic response.

% Water intake
fix_win     = false; % Boolean for fixing water intake.
low_win_ind = false; % Boolean for having low water intake.

% Sex specific mechanisms
m_RSNA = false; % Boolean for having male RSNA in the female model.
m_AT2R = false; % Boolean for having male effect of AT2R in the female model.

%% Read and assign optional parameters.

% The odd inputs of varargin are strings for each scenario. The
% corresponding even inputs are the values for the effect parameters to
% modify something.

kappa_ACEi_inp = 0;
kappa_ARB_inp = 0;

for i = 1:2:length(varargin)
    % Drugs
    if strcmp(varargin{i},'ACEi')
        kappa_ACEi_inp  = varargin{i + 1};
    elseif strcmp(varargin{i},'ARB')
        kappa_ARB_inp   = varargin{i + 1};
    end
end
kappa_ACEi  = kappa_ACEi_inp;
kappa_ARB = kappa_ARB_inp;

%% Retrieve parameters by name.

N_rsna           = pars(3 );
R_aass           = pars(4 );
R_eass           = pars(5 );
P_B              = pars(6 );
P_go             = pars(7 );
C_gcf            = pars(8 );
eta_ptsodreab_eq = pars(9 );
eta_dtsodreab_eq = pars(10);
eta_cdsodreab_eq = pars(11);
K_vd             = pars(12);
K_bar            = pars(13) * SF_R;
R_bv             = pars(14) * SF_R;
N_adhs_eq        = pars(15);
T_adh            = pars(16);
Phi_sodin_const  = pars(17);
N_als_eq         = pars(18);
C_K              = pars(19);
T_al             = pars(20);
N_rs             = pars(21);
X_PRCPRA         = pars(22);
h_renin          = pars(23);
h_AGT            = pars(24);
h_AngI           = pars(25);
h_AngII          = pars(26);
h_Ang17          = pars(27);
h_AngIV          = pars(28);
h_AT1R           = pars(29);
h_AT2R           = pars(30);
k_AGT            = pars(31);
c_ACE            = pars(32);
c_Chym           = pars(33);
c_NEP            = pars(34);
c_ACE2           = pars(35);
c_IIIV           = pars(36);
c_AT1R           = pars(37);
c_AT2R           = pars(38);
AT1R_eq          = pars(39);
AT2R_eq          = pars(40);
Psi_AT2RAA_eq    = pars(41);
Psi_AT2REA_eq    = pars(42);
ALD_eq           = pars(43);
% new calcium-magnesium parameters
C_Ca              = pars(44); % calcium concentration
C_Mg              = pars(45); % magnesium concentration
C_D3              = pars(46); % 1,25(OH)2D3 concentration
C_PTH             = pars(47); % PTH concentration
ratio_CaMg_eq     = pars(48);
gamma_r_CaMg      = pars(49);
K_sod_Ca          = pars(50);
K_sod_PTH         = pars(51);
K_ALD_Ca          = pars(52);
h_ALD_Ca          = pars(53);
K_ALD_Mg          = pars(54);
h_ALD_Mg          = pars(55);
K_ALD_PTH         = pars(56);
h_ALD_PTH         = pars(57);
K_R_Ca            = pars(58);
h_R_Ca            = pars(59);
K_R_D3_l            = pars(60);
K_R_D3_h            = pars(61);

% Species specific parameters for water transport
if     strcmp(species, 'human')
pt_sod_reab_EQ   = pars(62);
dt_sod_reab_EQ   = pars(63);
cd_sod_reab_EQ   = pars(64);
A_twreab         = pars(65);
elseif strcmp(species, 'rat')
eta_ptwreab_eq   = pars(62);
eta_dtwreab_eq   = pars(63);
eta_cdwreab_eq   = pars(64);
end

%% Retrieve variables by name.

rsna           = x(1 ); 
alpha_map      = x(2 ); 
alpha_rap      = x(3 );
alpha_Ca       = x(4 );
alpha_Mg       = x(5 );
R_r            = x(6 );
beta_rsna      = x(7 ); 
Phi_rb         = x(8 ); 
Phi_gfilt      = x(9 ); 
P_f            = x(10 );  
P_gh           = x(11); 
Sigma_tgf      = x(12); 
Phi_filsod     = x(13); 
Phi_ptsodreab  = x(14); 
eta_ptsodreab  = x(15); 
gamma_filsod   = x(16); 
gamma_at       = x(17); 
gamma_rsna     = x(18);
gamma_casr     = x(19);
gamma_PTH      = x(20);
Phi_mdsod      = x(21); 
Phi_dtsodreab  = x(22); 
eta_dtsodreab  = x(23); 
psi_al         = x(24);
Phi_dtsod      = x(25); 
Phi_cdsodreab  = x(26); 
eta_cdsodreab  = x(27); 
lambda_dt      = x(28); 
lambda_anp     = x(29); 
lambda_al      = x(30); 
Phi_usod       = x(31); 
Phi_sodin      = x(32); 
V_ecf          = x(33); 
V_b            = x(34); 
P_mf           = x(35); 
Phi_vr         = x(36); 
Phi_co         = x(37); 
P_ra           = x(38); P_ra_p = x_p(38); 
vas            = x(39); 
vas_f          = x(40); 
vas_d          = x(41);
Psi_AT1R_AngII = x(42);
Psi_Ca         = x(43);
Psi_Mg         = x(44);
Psi_PTH        = x(45);
Psi_D3         = x(46);
R_a            = x(47); 
R_ba           = x(48); 
R_vr           = x(49); 
R_tp           = x(50);  
P_ma           = x(51); 
epsilon_aum    = x(52); 
a_auto         = x(53); a_auto_p = x_p(53);
a_chemo        = x(54); 
a_baro         = x(55); 
C_adh          = x(56); 
N_adh          = x(57); 
N_adhs         = x(58); 
delta_ra       = x(59); 
M_sod          = x(60); 
C_sod          = x(61); 
nu_mdsod       = x(62); 
nu_rsna        = x(63); 
C_al           = x(64); 
N_al           = x(65); 
N_als          = x(66); 
xi_ksod        = x(67); 
xi_map         = x(68); 
xi_at          = x(69); 
xi_Ca          = x(70);
xi_Mg          = x(71);
xi_PTH         = x(72);
hatC_anp       = x(73);
AGT            = x(74);
nu_AT1         = x(75);
nu_Ca          = x(76);
nu_D3          = x(77);
R_sec          = x(78);
PRC            = x(79); 
PRA            = x(80); 
AngI           = x(81); 
AngII          = x(82); 
AT1R           = x(83);  
AT2R           = x(84);
Ang17          = x(85); 
AngIV          = x(86); 
R_aa           = x(87); 
R_ea           = x(88);
Sigma_myo      = x(89);
Psi_AT1RAA     = x(90); 
Psi_AT1REA     = x(91); 
Psi_AT2RAA     = x(92); 
Psi_AT2REA     = x(93);

% Species specific variables for water transport
if     strcmp(species, 'human')
Phi_twreab     = x(94);
mu_adh         = x(95); 
mu_Na          = x(96); 
Phi_u          = x(97); 
Phi_win        = x(98);
elseif strcmp(species, 'rat')
Phi_ptwreab    = x(94); 
eta_ptwreab    = x(95); 
mu_ptsodreab   = x(96); 
Phi_mdu        = x(97); 
Phi_dtwreab    = x(98); 
eta_dtwreab    = x(99);  
mu_dtsodreab   = x(100); 
Phi_dtu        = x(101);  
Phi_cdwreab    = x(102); 
eta_cdwreab    = x(103); 
mu_cdsodreab   = x(104); 
mu_adh         = x(105);  
Phi_u          = x(106); 
Phi_win        = x(107); 
end

% new calcium-magnesium variables

%% Differential algebraic equation system f(t,x,x') = 0.

f = zeros(length(x),1);

% rsna
rsna0 = N_rsna * alpha_map * alpha_rap * alpha_Ca * alpha_Mg;
if     strcmp(sex,'male')
        f(1 ) = rsna - rsna0;
elseif strcmp(sex,'female')
    if   m_RSNA
        f(1 ) = rsna - rsna0;
    else
        f(1 ) = rsna - rsna0^(1/rsna0);
    end
end

% alpha_map
if     strcmp(species, 'human')
    alphamap_a = 100;
elseif strcmp(species, 'rat')
    alphamap_a = fixed_var_pars(1);
end
f(2 ) = alpha_map - ( 0.5 + 1 / (1 + exp((P_ma - alphamap_a) / 15)) );

% alpha_rap
f(3 ) = alpha_rap - ( 1 - 0.008 * P_ra );

% alpha_Ca
f(4 ) = alpha_Ca - ( 0.85 + 0.3 / (1 + (1.232/C_Ca)) );

% alpha_Mg
f(5 ) = alpha_Mg - ( 0.85 + 0.3 / (1 + (C_Mg/0.64)) );

% R_r
f(6 ) = R_r - ( R_aa + R_ea );

% beta_rsna
if   denerve
    f(7 ) = beta_rsna - ( SSdata_fix(5) );
else
    f(7 ) = beta_rsna - ( 2 / (1 + exp(-3.16 * (rsna - 1))) );
end

% Phi_rb
if     RPP_ind
    f(8 ) = Phi_rb - ( RPP / R_r );
else
    f(8 ) = Phi_rb - ( P_ma / R_r );
end

% Phi_gfilt
f(9 ) = Phi_gfilt - ( P_f * C_gcf );

% P_f
f(10 ) = P_f - ( P_gh - (P_B + P_go) );

% P_gh
if     RPP_ind
    f(11) = P_gh - ( RPP - Phi_rb * R_aa );
else
    f(11) = P_gh - ( P_ma - Phi_rb * R_aa );
end

% Sigma_tgf
if     strcmp(species, 'human')
    Sigmatgf_a = 3.859 * SF_S;
elseif strcmp(species, 'rat')
    Sigmatgf_a = fixed_var_pars(2);
end
if     no_TGF
    f(12) = Sigma_tgf - ( 1 );
elseif NSAID == 1
    f(12) = Sigma_tgf - ( 0.644032 + 1.188073289 / (2.0285174154 + exp(((1-kappa_f_md) * Phi_mdsod - 3.859*SF_S)/(-0.9617*SF_S))));
elseif NSAID > 1
    f(12) = Sigma_tgf - ( 1 );
else
    f(12) = Sigma_tgf - ( 0.3408 + 3.449 / (3.88 + exp(((1-kappa_f_md) * Phi_mdsod - Sigmatgf_a) / (-0.9617 * SF_S) )) );
end

% Phi_filsod
f(13) = Phi_filsod - ( Phi_gfilt * C_sod );

% Phi_ptsodreab
f(14) = Phi_ptsodreab - ( Phi_filsod * eta_ptsodreab );

% eta_ptsodreab
f(15) = eta_ptsodreab - ( (1-kappa_f)* ((eta_ptsodreab_eq-0.20)*gamma_PTH + 0.20*gamma_casr) * gamma_filsod * gamma_at * gamma_rsna );

% gamma_filsod
if     strcmp(species, 'human')
    gammafilsod_a = 0.85; gammafilsod_b = 18;
elseif strcmp(species, 'rat')
    gammafilsod_a = 0.8 ; gammafilsod_b = fixed_var_pars(3);
end
f(16) = gamma_filsod - ( gammafilsod_a + 0.3 / (1 + exp((Phi_filsod - gammafilsod_b)/(138 * SF_S) )) );

% gamma_at
if     strcmp(species, 'human')
    gammaat_c = 2.6/2.3418; gammaat_d = 0.95; gammaat_a = 0.12 ; 
    gammaat_b = 2.3418;
elseif strcmp(species, 'rat')
    gammaat_c = 0.8017; gammaat_d = 0.92; gammaat_a = 0.136;
    gammaat_b = -1/(1-gammaat_c) * log(gammaat_a/(1-gammaat_d) - 1);
end
f(17) = gamma_at - ( gammaat_d + gammaat_a / (1 + exp(-gammaat_b * ((AT1R/AT1R_eq) - gammaat_c)) ) );

% gamma_rsna
if   denerve
    f(18) = gamma_rsna - ( SSdata_fix(16) );
else
    f(18) = gamma_rsna - ( 0.72 + 0.56 / (1 + exp((1 - rsna) / 2.18)) );
end

% gamma_CaSR
f(19) = gamma_casr - ( 0.904311039484287 + 0.19 / (1 + (C_Ca/K_sod_Ca)) );

% gamma_PTH
f(20) = gamma_PTH - ( 0.979389419757292 + 0.17 / (1 + (C_PTH/K_sod_PTH)^2) );

% Phi_mdsod
f(21) = Phi_mdsod - ( Phi_filsod - Phi_ptsodreab );

% Phi_dtsodreab
f(22) = Phi_dtsodreab - ( Phi_mdsod * eta_dtsodreab );

% eta_dtsodreab
f(23) = eta_dtsodreab - ( eta_dtsodreab_eq * psi_al * 1 );

% psi_al
if     strcmp(species, 'human')
    f(24) = psi_al - (0.17 + 0.94/(1+exp((0.48 - 1.2*log(C_al))/0.88)));
elseif strcmp(species, 'rat')
    psial_b = 0.1; psial_d = 1.05 / psial_b; psial_a = (1 + psial_b) * psial_d;
    psial_c = -1/ALD_eq * log((psial_a / (1 + psial_d) - 1) / psial_b);
    f(24) = psi_al - ( psial_a / (1 + psial_b * exp(-psial_c * C_al)) - psial_d );
end

% Phi_dtsod
f(25) = Phi_dtsod - ( Phi_mdsod - Phi_dtsodreab );

% Phi_cdsodreab
f(26) = Phi_cdsodreab - ( Phi_dtsod * eta_cdsodreab );

% eta_cdsodreab
eta_cdsodreab0 = ( eta_cdsodreab_eq * lambda_dt * lambda_anp * lambda_al);
if     strcmp(species, 'human')
    f(27) = eta_cdsodreab - ( eta_cdsodreab0 );
elseif strcmp(species, 'rat')
    etacd_a = 0.02;
    if     eta_cdsodreab0 <= eta_cdsodreab_eq + etacd_a
        f(27) = eta_cdsodreab - ( eta_cdsodreab0);
    elseif eta_cdsodreab0 > eta_cdsodreab_eq + etacd_a
        etacd_b = 12;
        etacd_c = (eta_cdsodreab_eq + etacd_a) - (1/etacd_b) * atanh((eta_cdsodreab_eq + etacd_a));
        f(27) = eta_cdsodreab - ( tanh(etacd_b * (eta_cdsodreab0 - etacd_c)) );
    end
end

% lambda_dt
if     strcmp(species, 'human')
    lambdadt_a = 0.82   ; lambdadt_b = 0.39         ;
    lambdadt_c = 1/0.375; lambdadt_d = 1.7625 * SF_S;
elseif strcmp(species, 'rat')
    lambdadt_a = 0.8;
    lambdadt_d = fixed_var_pars(4);
    lambdadt_b = 0.275; lambdadt_c = 2.3140;
end
f(28) = lambda_dt - ( lambdadt_a + lambdadt_b/ (1 + exp(lambdadt_c/SF_S * (Phi_dtsod - lambdadt_d)) ) );

% lambda_anp
f(29) = lambda_anp - ( -0.1 * hatC_anp + 1.1 );

% lambda_al
if     strcmp(species, 'human')
    lambdaal_a = 1/0.76;
elseif strcmp(species, 'rat')
    lambdaal_a = (ALD_eq^0.06);
end
f(30) = lambda_al - ( 1/lambdaal_a * C_al^0.06 );

% Phi_usod
if     strcmp(species, 'human')
    f(31) = Phi_usod - max(0,( Phi_dtsod - Phi_cdsodreab ));
elseif strcmp(species, 'rat')
    f(31) = Phi_usod -       ( Phi_dtsod - Phi_cdsodreab  );
end

% Phi_sodin
if     strcmp(species, 'human')
    sodin_E = 0.14; sodin_C = 0.14; sodin_D = 0.765; sodin_L = 0.1;
    sodin_B = (sodin_L*sodin_C - sodin_E*sodin_C +(sodin_E - 0.126)*sodin_C*ALD_eq^sodin_D)/(0.126-sodin_L);
    sodin_A = (0.126-sodin_E)*(sodin_B+sodin_C*ALD_eq^sodin_D);
    f(32) = Phi_sodin - max(0,sodin_A/(sodin_B+sodin_C*C_al^sodin_D)+sodin_E);
elseif strcmp(species, 'rat')
    f(32) = Phi_sodin - ( Phi_sodin_const );
end

% V_ecf
f(33) = Phi_win - Phi_u;

% V_b
if strcmp(species, 'human')
    if strcmp (sex, 'male')
        f(34) = V_b - ( 0.325 * V_ecf );
    else
        f(34) = V_b - ( 0.325 * V_ecf );
    end
elseif strcmp(species, 'rat')
    f(34) = V_b - ( SF_V*( 4.5479 + 2.4312 / (1 + exp(-(V_ecf - 18.1128*SF_V) * (0.4744/SF_V) )) ) );
end

% P_mf
if strcmp(sex,'male')
    pmf_a = (7.4360/SF_V);
else
    pmf_a = 1.1*(7.4360/SF_V);
end
f(35) = P_mf - ( ( pmf_a * V_b - 30.18) * epsilon_aum );

% Phi_vr
f(36) = Phi_vr - ( (P_mf - P_ra) / R_vr );

% Phi_co
f(37) = Phi_co - ( Phi_vr );

% P_ra
if     strcmp(species, 'human')
    if     strcmp(sex,'male')
        pra_a =  0; %0.8268;
    elseif strcmp(sex,'female')
        pra_a =   0; %0.8245;
    end
    f(38) = P_ra - ( 0.2787 * exp(Phi_co * 0.2281 * SF_R) - pra_a );
elseif strcmp(species, 'rat')
    pra_a = 0.2787 * exp(SSdata_input(37) * 0.2281 * SF_R);
    pra_a = pra_a + 0.01 * pra_a;
    f(38) = P_ra - ( max( 0, 0.2787 * exp(Phi_co * 0.2281 * SF_R) - pra_a ) );
end

% vas
f(39) = 1 / 1000 * (vas_f - vas_d);

% vas_f
if     strcmp(species, 'human')
    if strcmp(sex, 'male')
        vasf_a = 0.4799;
    else
        vasf_a = 0.4799; %*1.15;
    end
elseif strcmp(species, 'rat')
    vasf_a = -1/SSdata_input(37) * log(1/11.312);
end
f(40) = vas_f - ( (11.312 * exp(-Phi_co * vasf_a)) / 100 );

% vas_d
f(41) = vas_d - ( vas * K_vd );

% Psi_AT1R_AngII
xiat_a = 0.6; xiat_b = 2.4 - xiat_a; xiat_c = 2.5;
xiat_d = 1 + 1/xiat_c * log(xiat_b/(1-xiat_a) - 1);
f(42) = Psi_AT1R_AngII - ( xiat_a + xiat_b / (1 + exp(-xiat_c * ((AT1R/AT1R_eq) - xiat_d))) );

% Psi_Ca
f(43) = Psi_Ca - ( 0.865 + 0.27 / (1 + (1.232/C_Ca)^1) );

% Psi_Mg
f(44) = Psi_Mg - ( 0.85 + 0.3 / (1 + (C_Mg/0.64)^1) );

% Psi_PTH
PTH_eq = 6.4237;
f(45) = Psi_PTH - ( 0.8 + 0.4 / (1 + (PTH_eq/C_PTH)^1) );

% Psi_D3
f(46) = Psi_D3 - ( 0.85 + 0.3 / (1 + (152.68/C_D3)^4) );

% R_a
f(47) = R_a - ( R_ba * epsilon_aum * Psi_AT1R_AngII * Psi_Ca * Psi_Mg * Psi_PTH );

% R_ba
f(48) = R_ba - ( K_bar / vas );

% R_vr
f(49) = R_vr - ( (8 * R_bv + R_a) / 31 );

% R_tp
f(50) = R_tp - ( R_a + R_bv );

% P_ma
f(51) = P_ma - ( Phi_co * R_tp );

% epsilon_aum
if     strcmp(species, 'human')
    epsilonaum_a = 1;
elseif strcmp(species, 'rat')
    epsilonaum_a = 4/5;
end
f(52) = epsilon_aum - ( epsilonaum_a * (a_chemo + a_baro) );

% a_auto
if     strcmp(species, 'human')
    aauto_a = 0.011;
elseif strcmp(species, 'rat')
    aauto_a = fixed_var_pars(5); 
end
f(53) = a_auto - ( 3.0042 * exp(-aauto_a * P_ma) ); 

% a_chemo
f(54) = a_chemo - ( 1/4 * a_auto );

% a_baro
f(55) = 3/4 * (a_auto_p - 0.0000667 * (a_baro - 1));

% C_adh
f(56) = C_adh - ( 4 * N_adh );

% N_adh
f(57) = 1/T_adh * (N_adhs - N_adh);

% N_adhs
if     strcmp(species, 'human')
    Nadhs_a = 141;
elseif strcmp(species, 'rat')
    Nadhs_a = fixed_var_pars(6);
end
f(58) = N_adhs - ( N_adhs_eq * (max( 0, C_sod - Nadhs_a) + max( 0, epsilon_aum - 1 ) - delta_ra) / 3 );

% delta_ra
f(59) = 0.2 * P_ra_p - 0.0007 * delta_ra;

% M_sod
f(60) = Phi_sodin - Phi_usod;

% C_sod
f(61) = C_sod - ( M_sod / V_ecf );

% nu_mdsod
if     strcmp(species, 'human')
    if     strcmp(sex,'male')
        numdsod_a = 1.731 * SF_S;
    elseif strcmp(sex,'female')
        numdsod_a = 1.637 * SF_S;
    end
elseif strcmp(species, 'rat')
    numdsod_a = fixed_var_pars(7);
end
if NSAID > 0
    f(62) = nu_mdsod - ( 0.2262 + 83.4095/ (77.6196 + exp(((1-kappa_f_md) * Phi_mdsod - numdsod_a) / (0.6056 * SF_S) )) );
else
    f(62) = nu_mdsod - ( 0.2262 + 28.04  / (11.56   + exp(((1-kappa_f_md) * Phi_mdsod - numdsod_a) / (0.6056 * SF_S) )) );
end

% nu_rsna
if     strcmp(species, 'human')
    nursna_a = 0.8667;
elseif strcmp(species, 'rat')
    nursna_a = 0.8662;
end
if   denerve
    f(63) = nu_rsna - ( SSdata_fix(54) );
else
    f(63) = nu_rsna - ( 1.822 - 2.056 / (1.358 + exp(rsna - nursna_a)) );
end

% C_al
if     strcmp(species, 'human')
    f(64) = C_al - ( max( 1, N_al * ALD_eq ) );
elseif strcmp(species, 'rat')
    f(64) = C_al - (         N_al * ALD_eq   );
end

% N_al
f(65) = 1/T_al * (N_als - N_al);

% N_als
f(66) = N_als - ( N_als_eq * xi_ksod * xi_map * xi_at * xi_Ca * xi_Mg * xi_PTH );

% xi_ksod
if     strcmp(species, 'human')
    C_K_ref = 4.2;
    f(67) = xi_ksod - ( max( 0, (C_K / C_sod) / (C_K_ref/144/(6+1)) - 6 ) );
elseif strcmp(species, 'rat')
    f(67) = xi_ksod - ( 5 / ( 1 + exp(0.265 * (C_sod/C_K - fixed_var_pars(8))) ) ); 
end

% xi_map
if P_ma <= 100
    f(68) = xi_map - ( (1/exp(-0.0425 * 100)) * exp(-0.0425 * P_ma) );
else
    f(68) = xi_map - ( 1 );
end

% xi_at
if     strcmp(species, 'human')
    xiat_a = 0.47; xiat_b = 2.4; xiat_c = 1.5*1.301/0.8;
    xiat_d = 2.82/0.8/xiat_c;
elseif strcmp(species, 'rat')
    xiat_a = 0.2; xiat_b = 1.7 - xiat_a; xiat_c = 4.9;
    xiat_d = 1 + 1/xiat_c * log(xiat_b/(1-xiat_a) - 1);
end
f(69) = xi_at - ( xiat_a + xiat_b / (1 + exp(-xiat_c * ((AT1R/AT1R_eq) - xiat_d)) ) );

% xi_Ca
f(70) = xi_Ca - ( 0.649932812921984 + C_Ca^h_ALD_Ca / (K_ALD_Ca^h_ALD_Ca + C_Ca^h_ALD_Ca) );

% xi_Mg
f(71) = xi_Mg - ( 1.56641190964696 * K_ALD_Mg^h_ALD_Mg / (K_ALD_Mg^h_ALD_Mg + C_Mg^h_ALD_Mg) );

% xi_PTH
f(72) = xi_PTH - ( 0.997323497808287 + C_PTH^h_ALD_PTH / (K_ALD_PTH^h_ALD_PTH + C_PTH^h_ALD_PTH) );

% hatC_anp
f(73) = hatC_anp - ( 7.4052 - 6.554 / (1 + exp(P_ra - 3.762)) ); 

% AGT
f(74) = k_AGT - PRA - log(2)/h_AGT * AGT;

% nu_AT1
if     strcmp(species, 'human')
    f(75) = nu_AT1 - ( 10^(0.0102 - 0.95 * log10(AT1R / AT1R_eq)) );
elseif strcmp(species, 'rat')
    f(75) = nu_AT1 - ( (AT1R / AT1R_eq)^(-0.95) );
end

% nu_Ca
f(76) = nu_Ca - ( K_R_Ca^h_R_Ca / (K_R_Ca^h_R_Ca + C_Ca^h_R_Ca) );

% nu_D3
h = (C_D3/180)^50 / (1 + (C_D3/180)^50);
D3_impact_PRA_low = 0.99925739212222 + 0.5 * K_R_D3_l^4 / (C_D3^4 + K_R_D3_l^4);
D3_impact_PRA_high = 1 + C_D3^4 / (C_D3^4 + K_R_D3_h^4);

f(77) = nu_D3 - ( (1-h) * D3_impact_PRA_low + h * D3_impact_PRA_high );

% R_sec
f(78) = R_sec - ( N_rs * nu_mdsod * nu_rsna * nu_AT1 * nu_Ca * nu_D3 ); %nu_Ca * nu_D3

% PRC
f(79) = R_sec - log(2)/h_renin * PRC;
% PRA
f(80) = PRA - ( PRC * X_PRCPRA );
% AngI
f(81) = PRA - ((1-kappa_ACEi) * c_ACE + c_Chym + c_NEP) * AngI - log(2)/h_AngI * AngI;
% AngII
f(82) = kappa_AngII + ((1-kappa_ACEi) * c_ACE + c_Chym) * AngI - (c_ACE2 + c_IIIV + (1-kappa_ARB) * c_AT1R + c_AT2R) * AngII - log(2)/h_AngII * AngII;
% AT1R
f(83) = (1-kappa_ARB) * c_AT1R * AngII - log(2)/h_AT1R * AT1R;
% AT2R
f(84) = c_AT2R * AngII - log(2)/h_AT2R * AT2R;
% Ang17
f(85) = c_NEP * AngI + c_ACE2 * AngII - log(2)/h_Ang17 * Ang17;
% AngIV
f(86) = c_IIIV * AngII - log(2)/h_AngIV * AngIV;
% R_aa
f(87) = R_aa - ( R_aass * beta_rsna * Sigma_tgf * Sigma_myo * Psi_AT1RAA ...
                * Psi_AT2RAA * Psi_Ca * Psi_Mg * Psi_PTH ); % * Psi_D3 

% R_ea
f(88) = R_ea - ( R_eass * Psi_AT1REA * Psi_AT2REA );

% Sigma_myo
if     strcmp(species, 'human')
    if imp_Myo_ind
        sigmamyo_a = 1.2; sigmamyo_b = 0.3;
    else
        sigmamyo_a = 0.9; sigmamyo_b = 1.0;
    end 
    sigmamyo_d = 0.9;
    sigmamyo_c = 9;
    sigmamyo_e = 62;
elseif strcmp(species, 'rat')
    sigmamyo_a = 0.75; sigmamyo_b = 1.2; 
    sigmamyo_d = 0.6;
    sigmamyo_c = sigmamyo_b / (1-sigmamyo_a) - 1;
    sigmamyo_e = fixed_var_pars(9);
end
if     no_Myo
    f(89) = Sigma_myo - ( 1 );
elseif lin_Myo
    f(89) = Sigma_myo - ( 5 * (P_gh / sigmamyo_e - 1) + 1 );
else
    f(89) = Sigma_myo - ( sigmamyo_a + sigmamyo_b / ( 1 + sigmamyo_c * exp(-sigmamyo_d * (P_gh - sigmamyo_e)) ) );
end

% Psi_AT1RAA
f(90) = Psi_AT1RAA - ( 0.8   + 0.2092 * (AT1R / AT1R_eq) - 0.0092 / (AT1R / AT1R_eq) );

% Psi_AT1REA
f(91) = Psi_AT1REA - ( 0.925 + 0.0835 * (AT1R / AT1R_eq) - 0.0085 / (AT1R / AT1R_eq) );

% Psi_AT2RAA
if     strcmp(species, 'human')
    psiat2raa_a = 0.9;
    psiat2raa_b = 1 - psiat2raa_a;
    psiat2raa_c = 1;
elseif strcmp(species, 'rat')
    psiat2raa_a = 0.75;
    psiat2raa_b = 1 - psiat2raa_a;
    psiat2raa_c = 0.15;
end
if     strcmp(sex,'male')
        f(92) = Psi_AT2RAA - ( 1 );
elseif strcmp(sex,'female')
    if   m_AT2R
        f(92) = Psi_AT2RAA - ( 1 );
    else
        f(92) = Psi_AT2RAA - ( Psi_AT2RAA_eq * (psiat2raa_a + psiat2raa_b * exp(-psiat2raa_c * (AT2R/AT2R_eq - 1))) );
    end
end
% Psi_AT2REA
if     strcmp(species, 'human')
    psiat2rea_a = 0.9;
    psiat2rea_b = 1 - psiat2rea_a;
    psiat2rea_c = 1;
elseif strcmp(species, 'rat')
    psiat2rea_a = 0.8;
    psiat2rea_b = 1 - psiat2rea_a;
    psiat2rea_c = 0.15;
end
if     strcmp(sex,'male')
        f(93) = Psi_AT2REA - ( 1 );
elseif strcmp(sex,'female')
    if   m_AT2R
        f(93) = Psi_AT2REA - ( 1 );
    else
        f(93) = Psi_AT2REA - ( Psi_AT2REA_eq * (psiat2rea_a + psiat2rea_b * exp(-psiat2rea_c * (AT2R/AT2R_eq - 1))) );
    end    
end


if     strcmp(species, 'human')
    % Phi_twreab
    if strcmp(sex, 'male')
        f(94) = Phi_twreab - ( A_twreab - 0.001* 1/(mu_adh) + (0.8 + 0.08 * tanh( 8.5 * (mu_Na-1)) ) * Phi_gfilt );
    else
        f(94) = Phi_twreab - ( 1.01*(A_twreab - 0.001* 1/(mu_adh)) + (0.8 + 0.08 * tanh( 8.5 * (mu_Na-1)) ) * Phi_gfilt );
    end
    % mu_adh
    f(95) = mu_adh - ( 0.3313 + 0.8 / (1 + exp( 0.6 - 3.7 * log10(C_adh) )) );
    % mu_Na
    f(96) = mu_Na -  ( (Phi_ptsodreab + Phi_dtsodreab + Phi_cdsodreab) / (pt_sod_reab_EQ + dt_sod_reab_EQ + cd_sod_reab_EQ) );
    % Phi_u
    if (kappa_ACEi == 0) && (kappa_f == 0) && (kappa_f_md == 0) && (NSAID == 0)
        f(97) = Phi_u - ( max(0.0003     , Phi_gfilt - Phi_twreab ) );
    else   
        f(97) = Phi_u - ( max(0.0003*0.2 , Phi_gfilt - Phi_twreab ) );
    end
    % Phi_win
    if low_win_ind
        f(98) = Phi_win - (max(0, 0.0177    / (3.9271  + 18.22*C_adh^-1.607) - 0.002));
    else
        if strcmp(sex,'male')
            f(98) = Phi_win - (max(0, 0.0078541 / (0.65451 + 18.22*C_adh^-1.607) - 0.002));
        else
            f(98) = Phi_win - (max(0, 0.9*(0.0078541 / (0.65451 + 18.22*C_adh^-1.607) - 0.002)));
        end
    end

elseif strcmp(species, 'rat')
    % Phi_ptwreab
    f(94) = Phi_ptwreab - ( Phi_gfilt * eta_ptwreab );
    % eta_ptwreab
    f(95) = eta_ptwreab - ( eta_ptwreab_eq * mu_ptsodreab );
    % mu_ptsodreab
    musodreab_a = 0.12; musodreab_b = 10;
    f(96) = mu_ptsodreab - ( musodreab_a * tanh(musodreab_b * (eta_ptsodreab/eta_ptsodreab_eq - 1)) + 1 );
    % Phi_mdu
    f(97) = Phi_mdu - ( Phi_gfilt - Phi_ptwreab );
    % Phi_dtwreab
    f(98) = Phi_dtwreab - ( Phi_mdu * eta_dtwreab );
    % eta_dtwreab
    f(99) = eta_dtwreab - ( eta_dtwreab_eq * mu_dtsodreab );
    % mu_dtsodreab
    f(100) = mu_dtsodreab - ( musodreab_a * tanh(musodreab_b * (eta_dtsodreab/eta_dtsodreab_eq - 1)) + 1 );
    % Phi_dtu
    f(101) = Phi_dtu - ( Phi_mdu - Phi_dtwreab );
    % Phi_cdwreab
    f(102) = Phi_cdwreab - ( Phi_dtu * eta_cdwreab );
    % eta_cdwreab
    f(103) = eta_cdwreab - ( eta_cdwreab_eq * mu_cdsodreab * mu_adh );
    % mu_cdsodreab
    f(104) = mu_cdsodreab - ( musodreab_a * tanh(musodreab_b * (eta_cdsodreab/eta_cdsodreab_eq - 1)) + 1 );
    % mu_adh
    muadh_a = 1.0328; muadh_b = 0.1938;
    muadh_c = -1/4 * log((muadh_a - 1) / muadh_b);
    f(105) = mu_adh - ( muadh_a - muadh_b * exp(-muadh_c * C_adh) );
    % Phi_u
    f(106) = Phi_u - ( Phi_dtu - Phi_cdwreab );
    
    % Phi_win
    if     fix_win
        f(107) = Phi_win - ( SSdata_fix(93) );
    else
        phiwin_a = 0.94; phiwin_c = 0.0027; phiwin_d = 0.025;
        phiwin_b = SSdata_input(56) + 1 / phiwin_a * log(phiwin_c*SF_U / (0.030 - phiwin_d) - 1);
        f(107) = Phi_win - ( (phiwin_c * SF_U / (1 + exp(-phiwin_a * (C_adh - phiwin_b))) + phiwin_d) );
    end
end

end
