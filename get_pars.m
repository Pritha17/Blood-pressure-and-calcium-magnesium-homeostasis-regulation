function pars = get_pars(species, sex, varargin)

%% Species and sex identifier

% Species boolean
if     strcmp(species, 'human')
    spe_par = 1;
elseif strcmp(species, 'rat')
    spe_par = 0;
end
% Sex boolean
if     strcmp(sex, 'male')
    sex_par = 1;
elseif strcmp(sex, 'female')
    sex_par = 0;
end


%% Default parameter inputs for changes in simulation.

m_Reab = false; % Boolean for having male fractional sodium/water reabsorption in the female model.
m_RAS  = false; % Boolean for having male RAS parameters in the female model.
f_Reab = false; % Boolean for having female fractional sodium/water reabsorption in the male model.
f_RAS  = false; % Boolean for having female RAS parameters in the male model.
htn_rsna= 1    ; % Multiplicative factor for RSNA to induce hypertension.
htn_renin  = 1    ; % Multiplicative factor for Renin secretion to induce hypertension.
htn_raa    = 1    ; % Multiplicative factor for baseline AA resistance to induce hypertension.
htn_ald    = 1    ;
low_gfr    = 1    ; % <1 to similar low GFR

%% Read and assign optional parameters.

% The odd inputs of varargin are strings for each scenario. The
% corresponding even inputs are the values for the effect parameters to
% modify something.

for i = 1:2:length(varargin)
    if     strcmp(varargin{i},'Normal')
    elseif strcmp(varargin{i},'m_Reab') || strcmp(varargin{i},'m_RSNA_m_Reab')
        m_Reab = varargin{i + 1};
    elseif strcmp(varargin{i},'f_Reab')
        f_Reab = varargin{i + 1};
    elseif strcmp(varargin{i},'m_RAS' )
        m_RAS  = varargin{i + 1};
    elseif strcmp(varargin{i},'f_RAS' )
        f_RAS  = varargin{i + 1};
    elseif strcmp(varargin{i},'m_RAS_m_Reab' )
        m_Reab = varargin{i + 1};
        m_RAS  = varargin{i + 1};
    elseif strcmp(varargin{i},'RSNA' )
        htn_rsna   = varargin{i + 1};
    elseif strcmp(varargin{i},'Renin' )
        htn_renin  = varargin{i + 1};
    elseif strcmp(varargin{i},'R_aa' )
        htn_raa    = varargin{i + 1};
    elseif strcmp(varargin{i},'ALD' )
        htn_ald    = varargin{i + 1};
    elseif strcmp(varargin{i},'low GFR' )
        low_gfr    = varargin{i + 1};
    end
end

%% Generic parameters

N_rsna     = 1 * htn_rsna; % -
P_B        = 18     ; % mmHg
P_go       = 28     ; % mmHg
K_vd       = 0.01   ; % -
if strcmp (sex, 'male') && not(f_RAS) && not(f_Reab)
    K_bar      = 16.6   ; % mmHg min / ml
    R_bv       = 3.4    ; % mmHg min / ml
else
    K_bar      = 16.6 *1.2   ; % mmHg min / ml
    R_bv       = 3.4  *1.2    ; % mmHg min / ml
end
N_adhs_eq  = 1      ; % -
T_adh      = 6      ; % min
N_als_eq   = 1 * htn_ald     ; % -
C_K        = 5    ; % microEq / ml 
T_al       = 60     ; % min
N_rs       = 1 * htn_renin     ; % ng / ml / min

% RAS
h_renin       = 12  ; % min
h_AGT         = 600 ; % min
h_AngI        = 0.5 ; % min
h_AngII       = 0.66; % min
h_Ang17       = 30  ; % min
h_AngIV       = 0.5 ; % min
h_AT1R        = 12  ; % min
h_AT2R        = 12  ; % min
Psi_AT2RAA_eq = 1   ; % -
Psi_AT2REA_eq = 1   ; % -

% new calcium-magnesium parameters
C_Ca          = 1.232012; % mmol/L, calcium concentration
C_Mg          = 0.6403819; % mmol/L, magnesium concentration
C_D3          = 152.678631; % pmol/L, 1,25(OH)2D3 concentration
C_PTH         = 6.4237327; % pmol/L, PTH concentration
C_Ca_eq       = 1.232012; % mmol/L, calcium concentration
C_Mg_eq       = 0.6403819; % mmol/L, magnesium concentration
ratio_CaMg_eq = C_Ca_eq / C_Mg_eq;
gamma_r_CaMg  = 0.2;
K_sod_Ca = 1.25; % mM
K_sod_PTH = 2.386; % pM
K_ALD_Ca = 2.67; % mM
h_ALD_Ca = 0.8;
K_ALD_Mg = 1.0284; % mM
h_ALD_Mg = 1.2;
K_ALD_PTH = 124; % pM
h_ALD_PTH = 2;
K_R_Ca = 5.4; % mM
h_R_Ca = 6;
K_R_D3_l = 30; % pM
K_R_D3_h = 200; % pM

%% Species and sex specific parameters

if     strcmp(species, 'human')
    if strcmp(sex, 'male')
        R_aass    = 31.67 * htn_raa * (1/low_gfr) ; % mmHg min / l
        R_eass    = 51.66  * (1/low_gfr) ; % mmHg min / l
        C_gcf     = low_gfr * 0.00781;
        Phi_sodin = 0.126  ;
    else
        R_aass    = 31.67 * htn_raa ; % mmHg min / l
        R_eass    = 51.66  ; % mmHg min / l
        C_gcf     = low_gfr * 0.00781 ; %*0.8;
        Phi_sodin = 0.126  ;
    end
    
    %Sodium transport
    % adjust based on GFR
    scale_reabs = 0.1+.9*low_gfr;
    if (strcmp(sex,'male') && not(f_Reab)) || m_Reab
        eta_ptsodreab_eq = scale_reabs * 0.80;
        eta_dtsodreab_eq = scale_reabs * 0.50; 
        eta_cdsodreab_eq = scale_reabs * 0.93;
    else
        eta_ptsodreab_eq = scale_reabs * 0.65;
        eta_dtsodreab_eq = scale_reabs * 0.63; 
        eta_cdsodreab_eq = scale_reabs * 0.93;
        % female with male-like DT contribution
%         eta_ptsodreab_eq = scale_reabs * 0.65;
%         eta_dtsodreab_eq = scale_reabs * 0.50; 
%         eta_cdsodreab_eq = scale_reabs * 0.948;
%         eta_ptsodreab_eq = scale_reabs * 0.80;
%         eta_dtsodreab_eq = scale_reabs * 0.50; 
%         eta_cdsodreab_eq = scale_reabs * 0.93;
    end

    %RAAS
    X_PRCPRA    = 61/60.0  ; %fmol/min/pg
    if  (strcmp(sex, 'male') && not(f_RAS)) || m_RAS
        k_AGT   = 577.04   ;
        c_ACE   = 0.88492  ;
        c_Chym  = 0.09315  ;
        c_NEP   = 0.038189 ;
        c_ACE2  = 0.0078009;
        c_IIIV  = 0.25056  ;
        c_AT1R  = 0.17008  ;
        c_AT2R  = 0.065667 ;
        AT1R_eq = 13.99    ;
        AT2R_eq = 5.0854   ;
        ALD_eq  = 85       ;
        if   htn_rsna * htn_renin * htn_raa > 1
            A_twreab = 0.0182;
        else
            A_twreab = 0.0193;
        end
    elseif strcmp(sex, 'female') || f_RAS
        k_AGT   = 610.39   ;
        c_ACE   = 1.4079   ;
        c_Chym  = 0.1482   ;
        c_NEP   = 0.060759 ;
        c_ACE2  = 0.0037603;
        c_IIIV  = 0.038644 ;
        c_AT1R  = 0.027089 ;
        c_AT2R  = 0.038699 ;
        AT1R_eq = 3.78     ;
        AT2R_eq = 5.0854   ;
        ALD_eq  = 69.1775  ;
        if   htn_rsna * htn_renin * htn_raa  > 1
            A_twreab = 0.0181;
        else
            A_twreab = 0.0199;
        end
    end
    
    pt_sod_reab_EQ = 13.909;
    dt_sod_reab_EQ = 1.5859;
    cd_sod_reab_EQ = 1.6909;
    
elseif strcmp(species, 'rat')
    Phi_sodin = 2.4424; % microEq / min
    ALD_eq    = 387   ; % ng / l
    
    if     strcmp(sex, 'male')
        R_aass = 2.9905 * htn_raa; % mmHg min / ml
        R_eass = 4.878; % mmHg min / ml
        C_gcf  = 0.136; % ml / min / mmHg
        
        % Transport parameters
        eta_ptsodreab_eq = 0.80; % %
        eta_dtsodreab_eq = 0.50; % %
        eta_cdsodreab_eq = 0.93; % %
        eta_ptwreab_eq   = 0.86; % %
        eta_dtwreab_eq   = 0.60; % %
        eta_cdwreab_eq   = 0.78; % %
        
        % RAS
        X_PRCPRA = 135.59/17.312   ; % 1 / min
        k_AGT    = 801.02          ; % fmol / ml / min
        c_ACE    = 0.096833        ; % 1 / min
        c_Chym   = 0.010833        ; % 1 / min
        c_NEP    = 0.012667        ; % 1 / min
        c_ACE2   = 0.0026667       ; % 1 / min
        c_IIIV   = 0.29800         ; % 1 / min
        c_AT1R   = 0.19700         ; % 1 / min
        c_AT2R   = 0.065667        ; % 1 / min
        AT1R_eq  = 20.4807902818665; % fmol / ml
        AT2R_eq  = 6.82696474842298; % fmol / ml
    elseif strcmp(sex, 'female')
        R_aass = 4.6805 * htn_raa; % mmHg min / ml
        R_eass = 7.635; % mmHg min / ml
        C_gcf  = 0.094; % ml / min / mmHg
        
        % Transport parameters
        if   m_Reab
        eta_ptsodreab_eq = 0.71; % % % male
        eta_dtsodreab_eq = 0.50; % %
        eta_cdsodreab_eq = 0.93; % %
        eta_ptwreab_eq   = 0.80; % % % male 
        eta_dtwreab_eq   = 0.60; % %
        eta_cdwreab_eq   = 0.78; % %
        else
        eta_ptsodreab_eq = 0.50; % %
        eta_dtsodreab_eq = 0.50; % %
        eta_cdsodreab_eq = 0.96; % %
        eta_ptwreab_eq   = 0.50; % %
        eta_dtwreab_eq   = 0.60; % %
        eta_cdwreab_eq   = 0.91; % %
        end
        
        % RAS
        if   m_RAS
        X_PRCPRA = 135.59/17.312   ; % 1 / min % male
        k_AGT    = 801.02          ; % fmol / ml / min
        c_ACE    = 0.096833        ; % 1 / min
        c_Chym   = 0.010833        ; % 1 / min
        c_NEP    = 0.012667        ; % 1 / min
        c_ACE2   = 0.0026667       ; % 1 / min
        c_IIIV   = 0.29800         ; % 1 / min
        c_AT1R   = 0.19700         ; % 1 / min
        c_AT2R   = 0.065667        ; % 1 / min
        AT1R_eq  = 20.4807902818665; % fmol / ml
        AT2R_eq  = 6.82696474842298; % fmol / ml
        else
        X_PRCPRA = 114.22/17.312   ; % 1 / min
        k_AGT    = 779.63          ; % fmol / ml / min
        c_ACE    = 0.11600         ; % 1 / min
        c_Chym   = 0.012833        ; % 1 / min
        c_NEP    = 0.0076667       ; % 1 / min
        c_ACE2   = 0.00043333      ; % 1 / min
        c_IIIV   = 0.29800         ; % 1 / min
        c_AT1R   = 0.19700         ; % 1 / min
        c_AT2R   = 0.065667        ; % 1 / min
        AT1R_eq  = 20.4538920068419; % fmol / ml
        AT2R_eq  = 6.81799861123497; % fmol / ml
        end
    end
end

%% Human to rat scaling factors

% Scaling factors
% Rat value = Human value x SF
% Note: This includes conversion of units.

% Physiological variables which determine scaling factors.
% Original values are added as separate parameters because these may be
% modified by another script.

if     strcmp(species, 'human')
    SF_S =                     1; % sodium flow
    SF_U =                     1; % urine flow
    SF_R =                     1; % resistance
    SF_V =                     1; % volume
elseif strcmp(species, 'rat')
    if     strcmp(sex, 'male')
        Phi_sodin_orig = 2.4424 ; % Munger 1988, Karaaslan 2005
        Phi_u_orig = 0.030      ; % Munger 1988, Layton 2016
        R_r_orig = 2.9905 + 4.878; % Munger 1988
        W_b = 238               ; % Munger 1988
    elseif strcmp(sex, 'female')
        Phi_sodin_orig = 2.4424 ; % Munger 1988, Karaaslan 2005
        Phi_u_orig = 0.030      ; % Munger 1988, Layton 2016
        R_r_orig = 4.6805 + 7.635; % Munger 1988
        W_b = 194               ; % Munger 1988
    end
    V_b_orig = 0.06 * W_b + 0.77; % Lee 1985
    
    SF_S = Phi_sodin_orig / 0.126; % sodium flow
    SF_U = Phi_u_orig    / 0.001; % urine flow
    SF_R = R_r_orig      / 83.3 ; % resistance
    SF_V = V_b_orig      / 5    ; % volume
end

%% Retrieve parameters in fixed variable equations (rat only).

% Add directory containing data.
mypath = pwd;
mypath = strcat(mypath, '/Rat_Data');
addpath(genpath(mypath))

% These are the shift parameters which ensure that effect variables are 1.
% Retrieve baseline steady state variable values.
% This is done if the species is rat.
% The female model with male reabsorption requires its own file to be
% loaded since changing the transport parameters will change the steady
% state value of several variables.
if strcmp(species, 'rat')
    % Set name for data file to be loaded based upon sex and scenario.
    if     m_Reab
        load_data_name1 = sprintf('%s_fixed_var_pars_scenario_m_Reab.mat', sex);
        load_data_name2 = sprintf('%s_%s_ss_data_scenario_m_Reab.mat'    , species,sex);
    else
        load_data_name1 = sprintf('%s_fixed_var_pars_scenario_Normal.mat', sex);
        load_data_name2 = sprintf('rat_male_ss_data_scenario_normal_combined.mat');
    end
    load(load_data_name1, 'fixed_var_pars');
    load(load_data_name2, 'SSdata');
end
SSdata_BP = SSdata(1:107);

%% Parameter output/input

% Generic model parameters
pars = [spe_par; sex_par; ...
        N_rsna; R_aass; R_eass; P_B; P_go; C_gcf; ...
        eta_ptsodreab_eq; eta_dtsodreab_eq; eta_cdsodreab_eq; ...     
        K_vd; K_bar; R_bv; ...
        N_adhs_eq; T_adh; Phi_sodin; N_als_eq; C_K; T_al; N_rs; ...
        X_PRCPRA; h_renin; h_AGT; h_AngI; h_AngII; h_Ang17; h_AngIV; ...
        h_AT1R; h_AT2R; k_AGT; c_ACE; c_Chym; c_NEP; c_ACE2; c_IIIV; ...
        c_AT1R; c_AT2R; AT1R_eq; AT2R_eq; ...
        Psi_AT2RAA_eq; Psi_AT2REA_eq; ALD_eq];

% new calcium-magnesium parameters
pars = [pars; C_Ca; C_Mg; C_D3; C_PTH; ratio_CaMg_eq; gamma_r_CaMg;... 
        K_sod_Ca; K_sod_PTH; K_ALD_Ca; h_ALD_Ca; K_ALD_Mg; h_ALD_Mg;...
        K_ALD_PTH; h_ALD_PTH; K_R_Ca; h_R_Ca; K_R_D3_l; K_R_D3_h]; % size = 61

% Species model specific parameters
if     strcmp(species, 'human')
    pars = [pars; ...
            pt_sod_reab_EQ; dt_sod_reab_EQ; cd_sod_reab_EQ; A_twreab; ...
            SF_S; SF_U; SF_R; SF_V];
elseif strcmp(species, 'rat')
    pars = [pars; ...
            eta_ptwreab_eq; eta_dtwreab_eq; eta_cdwreab_eq; ...
            SF_S; SF_U; SF_R; SF_V; ...
            fixed_var_pars; SSdata_BP];
end

end
