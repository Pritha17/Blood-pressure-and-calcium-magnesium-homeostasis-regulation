% computes local sensitivity analysis of all parameters

incr = 1; % increase (1) or decrease (0) parameters by 5%

SSfile = './Rat_Data/rat_male_ss_data_scenario_normal_combined.mat';

save_res = 1;

dat = load(SSfile);
SS_IG = dat.SSdata;

pars_name_BP = ["R_{aa-eq}"; "R_{ea-eq}"; ... 
                "\eta_{pt+tal-sodreab}^{eq}"; "\eta_{dt-sodreab}^{eq}"; "\eta_{cd-sodreab}^{eq}";... 
                "[AT1R]_{eq}"; "[AT2R]_{eq}"; "[ALD]_{eq}"; "\gamma_r^{CaMg}"; "K_{sod-Ca}";...
                "K_{sod-PTH}"; "K_{ALD-Ca}"; "h_{ALD-Ca}"; "K_{ALD-Mg}"; "h_{ALD-Mg}";...
                "K_{ALD-PTH}"; "h_{ALD-PTH}"; "K_{R-Ca}"; "h_{R-Ca}"; "K_{R-D3-low}"; "K_{R-D3-high}"; ...
                "\eta_{pt-wreab}^{eq}"; "\eta_{dt-wreab}^{eq}"; "\eta_{cd-wreab}^{eq}"]; % 24

pars_name_CaMg = ["\gamma_p"; "\beta_{exo}^{PTHg}"; "k_{prod}^{PTHg}";... 
            "k_{deg}^{PTHg}"; "k_{deg}^{PTHp}"; "D_{3}^{inact}"; "\gamma_{deg}^{PTHp}"; "k_{deg}^{D3}";... 
            "I_{Ca}"; "\delta_{abs}^{D3}"; "I_{Mg}"; "V_{active}"; "K_{active}";... 
            "\lambda_{Ca-PT}^{0}"; "\lambda_{Ca-TAL}^{0}"; "\lambda_{Ca-DT}^{0}";... 
            "\lambda_{Mg-PT}^{0}"; "\lambda_{Mg-TAL}^{0}"; "\lambda_{Mg-DT}^{0}";... 
            "\delta_{res}^{max}"; "k_{pf}^{Ca}"; "k_{fp}^{Ca}"; "\gamma_{ac}^{Ca}";... 
            "k_{pf}^{Mg}"; "k_{fp}^{Mg}"; "\gamma_{ac}^{Mg}"; "K_{PTH-ALD}"; "K_{PTH-AT1R}"]; % 28

pars_name = [pars_name_BP; pars_name_CaMg];

species = 'rat';
sex = 'male';
[pars_BP, pars_CaMg, M] =  get_params_and_mass_matrix(SS_IG,species,sex,1,1,1,1);

pars = [pars_BP; pars_CaMg];

SSbase = compute_SS(pars, SS_IG, -1, incr, M);

indices = [4:5 9:11 39:40 43 49:64 185:186 189 193:195 197:198 208 210 ...
          212 216:217 219 221 223 225 227 229 240 245:250 265:266];
sens = zeros(length(indices), 9);
indx = [1 9 51 64 80 109:112];
frac_change = zeros(size(sens));
basevals = SSbase(indx)';

for ii = 1:length(indices)
    SSsens = compute_SS(pars, SS_IG, indices(ii), incr, M);
    sens(ii, :) = SSsens(indx);
    
    for jj = 1:length(indx)
        if basevals(jj) <= 1e-10
            frac_change(ii,jj) = 0;
        else
            frac_change(ii,jj) = 100.0*(sens(ii,jj) - basevals(jj))./basevals(jj);
        end
    end
end

%% save results
if save_res
    if incr
        type_change = 'I';
    else
        type_change = 'D';
    end
    fname = strcat('./Rat_Data/localsensitivity_male_',type_change,'.mat');
    save(fname, 'pars', 'frac_change', 'SSbase', 'sens', 'pars_name');
    fprintf('sensitivity analysis results saved to %s \n', fname)
end

function SS = compute_SS(pars, SS0, parchange_ID, incr, M)
    params = pars;
    if parchange_ID > -1
        if incr == 1
            f_change = 1.05;
        else
            f_change = 0.95;
        end
        params(parchange_ID) = f_change * pars(parchange_ID);
    end
    
    params_BP = params(1:184);
    params_CaMg = params(185:266);
    
    tchange=0;
    tspan = [0 1000000];
    options = odeset('Mass',M,'RelTol',1e-2,'AbsTol',1e-3*ones(1,length(SS0)));

    [t,x] = ode15s(@(t,x) all_eqns_bp_Mg(t,x,params_BP,params_CaMg,tchange,... 
                                            'ACEi',0, 'ARB', 0),...
                                            tspan,SS0, options);
                                          
    y=x(end,:);
    SS = y';
end
