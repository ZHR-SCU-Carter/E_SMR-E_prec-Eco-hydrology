
data_ALT_orig = data_ALT_orig;
data_prec_orig = data_prec_orig;
data_SM_root_orig = data_SM_root_orig;
data_SM_bottom_orig = data_SM_bottom_orig;

raw_data = [data_ALT_orig, data_prec_orig, data_SM_root_orig, data_SM_bottom_orig];
raw_data(any(isnan(raw_data), 2), :) = []; 

ALT  = raw_data(:, 1);
Prec = raw_data(:, 2);
SM_types = {"Root-zone SM", "Bottom-layer SM"};
SM_data  = [raw_data(:, 3), raw_data(:, 4)]; 

alt_threshold = 2.2;
groups = {ALT < alt_threshold, ALT >= alt_threshold};
group_names = {"ALT < 2.2m", "ALT >= 2.2m"};

fprintf('============================================================\n');

for g = 1:2 
    idx = groups{g};
    curr_ALT  = ALT(idx);
    curr_Prec = Prec(idx);
    
    X_ALT  = zscore(curr_ALT);
    X_Prec = zscore(curr_Prec);
    X_Full = [X_ALT, X_Prec]; 
    
    for s = 1:2 
        curr_SM = SM_data(idx, s);
        
        lm_full = fitlm(X_Full, curr_SM);
        R2_full = lm_full.Rsquared.Ordinary;
        
        p_model = coefTest(lm_full); 
        lm_noALT = fitlm(X_Prec, curr_SM);
        R2_noALT = lm_noALT.Rsquared.Ordinary;

        lm_noPrec = fitlm(X_ALT, curr_SM);
        R2_noPrec = lm_noPrec.Rsquared.Ordinary;

        unique_ALT  = max(0, R2_full - R2_noALT); 
        unique_Prec = max(0, R2_full - R2_noPrec);

        sum_unique = unique_ALT + unique_Prec;
        if sum_unique > 0
            rc_ALT  = (unique_ALT / sum_unique) * 100;
            rc_Prec = (unique_Prec / sum_unique) * 100;
        else
            rc_ALT = 0; rc_Prec = 0;
        end
        depth_name = SM_types{s};
        grp_name   = group_names{g};
        fprintf(' %-10s | Depth: %-15s (n=%d)\n', grp_name, depth_name, sum(idx));
        fprintf('  R² = %.4f | Full Model P = %.4e\n', R2_full, p_model);
        fprintf('  ALT Contri.: %6.2f%% | Prec Contri.: %6.2f%%\n', rc_ALT, rc_Prec);
        fprintf('------------------------------------------------------------\n');
    end
end
