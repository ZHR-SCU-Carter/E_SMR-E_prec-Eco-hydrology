%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a matlab code used to obtain the theta_cri of each site. 
% The required data are the average daily sensible heat flux, latent heat flux and the moisture content of the shallow surface soil. 
% The calculation formula can be found in equation (2) in the main text
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% EF calculation

LE = LE; % LE is latent heat flux
H = H; % H is sensible heat flux
EF = LE./(LE+H);

%% theta_cri determination
SM = SM;
validIdx = ~isnan(SM) & ~isnan(EF);  % Detect the index of non-nan data
SM = SM(validIdx);
EF = EF(validIdx);
% Sort the data in ascending order of x
[SM_sorted, sortIdSM] = sort(SM); 
EF_sorted = EF(sortIdSM); % Reorder         

fun = @(params, x) (x < params(4)) .* (params(1) * x + params(2)) + (x >= params(4)) .* (params(1) * params(4) + params(2));
params0 = [a0, b0, c0, d0]; % Initial value
lb = [a1, b1, c1, d1]; % Lower limit
ub = [a2, b2, c2, d2]; % Upper limit
% Fit using lsqcurvefit
options = optimset('MaxIter', 1000, 'MaxFunEvals', 2000);
[params_fit, resnorm] = lsqcurvefit(fun, params0, SM_sorted, EF_sorted, lb, ub, options);
x_highres = linspace(min(SM_sorted), max(SM_sorted), 1000); 
y_fit_highres = fun(params_fit, x_highres);              

% plot
figure('OuterPosition',[737,854.6,627.2,276.8]);
hold on;
scatter(SM_sorted, EF_sorted, 'b');
plot(x_highres, y_fit_highres, 'r-', 'LineWidth', 2);
xline(params_fit(4), '--k', 'LineWidth', 2);
xlabel('SM (m^3m^-^3)');
ylabel('EF');
hold off;
box on
