clc;
clear all

E_SMR = E_SMR;
E_prec = E_prec;
SOS = SOS;
airT = airT;
Elevation = Elevation;
sand_pct = sand_pct;
silt_pct = silt_pct;
clay_pct = clay_pct;
vars = {E_SMR, E_prec, SOS_raw, Elevation, airT, sand_pct, silt_pct, clay_pct};
for i=1:length(vars), v=vars{i}; v(v==-9999)=NaN; vars{i}=v; end
[E_SMR, E_prec, SOS_raw, Elevation, airT, sand_pct, silt_pct, clay_pct] = vars{:};

Nrows = length(E_SMR);

% Suppose there are three regions, observed for 10 years, 12 years and 9 years respectively
n1 = 10;
n2 = 12;
n3 = 9;
regionID = [ repelem("S1",n1).'; ...
           repelem("S2",n2).'; ...
           repelem("S3",n3).'; ...
           ];

regionCat = categorical(regionID);
E_SMR_z = zscore(E_SMR, 0, 'omitnan');
E_prec_z = zscore(E_prec, 0, 'omitnan');
Elevation_z = zscore(Elevation, 0, 'omitnan');
airT_z = zscore(airT, 0, 'omitnan');

tex = [sand_pct, silt_pct, clay_pct];
[coeff, score, latent, ~, explained] = pca(tex, 'Rows','complete'); 
texturePC1_z = zscore(score(:,1), 0, 'omitnan');
fprintf('    Texture PCA1 explains %.1f%% variance. PCA loadings (sand,silt,clay) = [%.3f, %.3f, %.3f]\n', ...
    explained(1), coeff(1,1), coeff(2,1), coeff(3,1));

tbl = table(SOS_centered, E_SMR_z, E_prec_z, Elevation_z, airT_z, texturePC1_z, siteCat, ...
    'VariableNames', {'SOS_c','E_SMR','E_prec','Elevation','airT','texturePC1','site'});
tbl.E_SMR2 = tbl.E_SMR.^2; tbl.E_SMR3 = tbl.E_SMR.^3;
tbl.E_prec2 = tbl.E_prec.^2; tbl.E_prec3 = tbl.E_prec.^3;
tbl.airT2 = tbl.airT.^2;   tbl.airT3 = tbl.airT.^3;

data_all = rmmissing(tbl);
formula_ext = ['SOS_c ~ 1 + ' ...
    'E_SMR + E_SMR2 + E_SMR3 + ' ...
    'E_prec + E_prec2 + E_prec3 + ' ...
    'airT + airT2 + airT3 + ' ...
    'E_SMR:E_prec + E_SMR2:E_prec + E_SMR:E_prec2 + ' ...
    'E_SMR:airT + E_SMR2:airT + E_SMR:airT2 + ' ...
    'E_prec:airT + E_prec2:airT + E_prec:airT2 + ' ...
    'Elevation + texturePC1 + (1|site)'];

fprintf('    (B) Fitting Extended Model...\n');
lme_ext = fitlme(data_all, formula_ext, 'Optimizer', 'quasinewton', 'FitMethod', 'REML');
aic_e = lme_ext.ModelCriterion.AIC;
bic_e = lme_ext.ModelCriterion.BIC;
[psi_e, mse_e] = covarianceParameters(lme_ext);
var_rand_e = psi_e{1}(1,1);
var_res_e = mse_e;
yhat_e = fitted(lme_ext, 'Conditional', false);
var_fix_e = var(yhat_e);
total_var_e = var_fix_e + var_rand_e + var_res_e;
r2_marg_e = var_fix_e / total_var_e;
r2_cond_e = (var_fix_e + var_rand_e) / total_var_e;
sigma_site_e = sqrt(var_rand_e);

fprintf('\n    ========================================================\n');
fprintf('    FULL MODEL COMPARISON SUMMARY\n');
fprintf('    ========================================================\n');
fprintf('    Metric                     | Baseline       | Extended\n');
fprintf('    --------------------------------------------------------\n');
fprintf('    AIC                        | %-14.2f | %-14.2f\n', aic_b, aic_e);
fprintf('    BIC                        | %-14.2f | %-14.2f\n', bic_b, bic_e);
fprintf('    R2_marg (Fixed)            | %-14.4f | %-14.4f\n', r2_marg_b, r2_marg_e);
fprintf('    R2_cond (Fixed+Random)     | %-14.4f | %-14.4f\n', r2_cond_b, r2_cond_e);
fprintf('    Random Effect SD (Site)    | %-14.4f | %-14.4f\n', sigma_site_b, sigma_site_e);
fprintf('    --------------------------------------------------------\n');
fprintf('    >>> Spatial Heterogeneity Explained: %.2f%%\n', (sigma_site_b - sigma_site_e) / sigma_site_b * 100);
fprintf('    ========================================================\n');

E = data_all.E_SMR; P = data_all.E_prec; T = data_all.airT;
Elev = data_all.Elevation; Tex = data_all.texturePC1;
Y = data_all.SOS_c; 

X_raw = [ ...
    E, P, T, ...                   % 1-3: Linear
    E.^2, E.*P, E.*T, ...          % 4-6: E2, EP, ET
    P.^2, P.*T, T.^2, ...          % 7-9: P2, PT, T2
    E.^3, E.^2.*P, E.^2.*T, ...    % 10-12: E3, E2P, E2T
    E.*P.^2, E.*T.^2, ...          % 13-14: EP2, ET2
    P.^3, P.^2.*T, P.*T.^2, ...    % 15-17: P3, P2T, PT2
    T.^3, ...                      % 18: T3
    Elev, Tex                      % 19-20: Controls
    ];

term_names = { ...
    'ESMR', 'Eprec', 'Tair', ...
    'ESMR2', 'ESMRxEprec', 'ESMRxTair', ...
    'Eprec2', 'EprecxTair', 'Tair2', ...
    'ESMR3', 'ESMR2xEprec', 'ESMR2xTair', ...
    'ESMRxEprec2', 'ESMRxTair2', ...
    'Eprec3', 'Eprec2xTair', 'EprecxTair2', ...
    'Tair3', ...
    'Elevation', 'texturePC1'};

R_raw = corr(X_raw);
VIF_before = diag(inv(R_raw));
Xz = zscore(X_raw, 0, 1); 
[Q, ~] = qr(Xz, 0);
r = size(Q,2);
R_Q = corr(Q);
VIF_after = diag(inv(R_Q));
tblQ = array2table([Y, Q], 'VariableNames', ['SOS', arrayfun(@(k) sprintf('Q%d',k),1:r,'Uni',0)]);
mdlQ = fitlm(tblQ, ['SOS ~ 1 + ' strjoin(tblQ.Properties.VariableNames(2:end),' + ')]);
b_hat = mdlQ.Coefficients.Estimate;
Covb = mdlQ.CoefficientCovariance;
R2_fixed = mdlQ.Rsquared.Ordinary; 
nb = 20000;
bsim = mvnrnd(b_hat, Covb, nb); 
bs_coeff = bsim(:, 2:end); 
contrib_Q = (bs_coeff.^2) ./ sum(bs_coeff.^2, 2) * R2_fixed;

contrib_mean = mean(contrib_raw, 1)';
contrib_sd = std(contrib_raw, 0, 1)'; 
fprintf('\n    ==============================================================================\n');
fprintf('    DETAILED FACTOR ANALYSIS (VIF & CONTRIBUTION with SD)\n');
fprintf('    ==============================================================================\n');
fprintf('    %-14s | VIF(Before) | VIF(After) | Contrib(%%) | SD(%%)\n', 'Term Name');
fprintf('    ------------------------------------------------------------------------------\n');
for i = 1:length(term_names)
    fprintf('    %-14s | %-11.2f | %-10.2f | %-10.2f | %-10.2f\n', ...
        term_names{i}, VIF_before(i), VIF_after(i), ...
        contrib_mean(i)*100, contrib_sd(i)*100);
end
fprintf('    ------------------------------------------------------------------------------\n');
fprintf('    * Note: "After" refers to data orthogonality. SD = Standard Deviation.\n');
fprintf('    ==============================================================================\n');

pred_names = {'E_{SMR}', 'E_{prec}', 'T_{air}', 'Elevation', 'Texture'};
means = mean(contrib_preds);
stds = std(contrib_preds, 0, 1);   % SD
ci_lo = prctile(contrib_preds, 2.5);
ci_hi = prctile(contrib_preds, 97.5);

fprintf('\n>>> Aggregated Contribution Results (Mean ± SD [95%% CI]):\n');
for k=1:5
    fprintf('    %s: %.2f%% ± %.2f%% [%.2f%% - %.2f%%]\n', ...
        pred_names{k}, means(k)*100, stds(k)*100, ci_lo(k)*100, ci_hi(k)*100);
end

fprintf('\n>>> [Step 4] Plotting...\n');
figure('Name', 'Joint Effects', 'Color', 'w', 'Position', [100, 100, 900, 700]);
airT_levels = [-1, 0, 1]; E_prec_levels = [-1, 0, 1];
colors = [0.2, 0.6, 0.8; 0.5, 0.5, 0.5; 0.8, 0.3, 0.3];
line_styles = {':', '--', '-'};

esmr_vals = linspace(min(data_all.E_SMR), max(data_all.E_SMR), 100)';
dummy_site = repmat(data_all.site(1), 100, 1);

hold on; idx=1; legend_entries={};
for i = 1:3
    for j = 1:3
        t_new = table();
        t_new.E_SMR = esmr_vals;
        t_new.E_prec = repmat(E_prec_levels(j), 100, 1);
        t_new.airT = repmat(airT_levels(i), 100, 1);
        t_new.Elevation = zeros(100,1); t_new.texturePC1 = zeros(100,1); t_new.site = dummy_site;
        
        t_new.E_SMR2 = t_new.E_SMR.^2; t_new.E_SMR3 = t_new.E_SMR.^3;
        t_new.E_prec2 = t_new.E_prec.^2; t_new.E_prec3 = t_new.E_prec.^3;
        t_new.airT2 = t_new.airT.^2; t_new.airT3 = t_new.airT.^3;
        
        y_pred_c = predict(lme_ext, t_new, 'Conditional', false);
        y_pred_real = y_pred_c + mean_SOS; 
        
        plot(esmr_vals, y_pred_real, 'Color', colors(i,:), 'LineStyle', line_styles{j}, 'LineWidth', 2);
        legend_entries{idx} = sprintf('T=%d, P=%d', airT_levels(i), E_prec_levels(j)); idx=idx+1;
    end
end
grid on; box on;
xlabel('Standardized E_{SMR}'); ylabel('Predicted SOS (DOY)');
title('(a) Joint Effects of Hydrothermal Factors');
legend(legend_entries, 'Location', 'eastoutside');
