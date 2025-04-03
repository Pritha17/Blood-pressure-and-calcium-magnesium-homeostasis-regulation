function [htn_rsna, htn_renin, htn_raa, htn_ald] = get_htn_factors(ht)
    % multiplicative factors for different types of HTN
    switch ht
        case 'All'
            htn_rsna   = 1.28  ; % Multiplicative factor for RSNA to induce hypertension.
            htn_renin  = 1.28  ; % Multiplicative factor for Renin secretion to induce hypertension.
            htn_raa    = 1.28  ; % Multiplicative factor for baseline AA resistance to induce hypertension.
            htn_ald    = 1.28  ;
        case 'RSNA'
            htn_rsna   = 1.54; %1.6  ; % Multiplicative factor for RSNA to induce hypertension.
            htn_renin  = 1.0    ; % Multiplicative factor for Renin secretion to induce hypertension.
            htn_raa    = 1.0; %1.1    ; % Multiplicative factor for baseline AA resistance to induce hypertension.
            htn_ald    = 1.0  ;
        case 'Renin'
            htn_rsna   = 1.0  ; % Multiplicative factor for RSNA to induce hypertension.
            htn_renin  = 6.0    ; % Multiplicative factor for Renin secretion to induce hypertension.
            htn_raa    = 1.0  ; % Multiplicative factor for baseline AA resistance to induce hypertension.else
            htn_ald    = 1.0  ;
        case 'RAA'
            htn_rsna   = 1.0  ; % Multiplicative factor for RSNA to induce hypertension.
            htn_renin  = 1.0  ; % Multiplicative factor for Renin secretion to induce hypertension.
            htn_raa    = 2.9  ; % Multiplicative factor for baseline AA resistance to induce hypertension.else
            htn_ald    = 1.0  ;
        case 'ALD'
            htn_rsna   = 1.0  ; % Multiplicative factor for RSNA to induce hypertension.
            htn_renin  = 1.0  ; % Multiplicative factor for Renin secretion to induce hypertension.
            htn_raa    = 1.0  ; % Multiplicative factor for baseline AA resistance to induce hypertension.else
            htn_ald    = 30.0  ;
        otherwise
            disp('Unknown HTN type')
    end
end