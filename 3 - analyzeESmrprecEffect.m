
function Results = analyzeESmrprecEffect(SM_check, time_line, years, len, ...
        theta_crisis, start_datee, end_datee, linking_date, theta_pf, pf_datetime, final_thawed_date)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This is the function for calculating E_SMR and E_prec. 
    % 
    % Based on the output of the function in "2 - extractAndVisualizeSMData.m" 
    % and manually specified inputs — theta_crisis, start_datee, end_datee, 
    % linking_date, theta_pf, pf_datetime, and final_thawed_date — the final 
    % results of E_SMR and E_prec are obtained.
    %
    % - start_datee: visually identified starting date for curve fitting.
    % 
    % - end_datee: visually identified ending date for curve fitting.
    % 
    % - linking_date: the date corresponding to the end point of the sliding 
    %   fitted curve.
    % 
    % - theta_pf: pre-freezing soil moisture for each freeze–thaw cycle.
    % 
    % - pf_datetime: the datetime timestamp associated with theta_pf.
    % 
    % - final_thawed_date: selected final thawing date, typically identified 
    %   as the point where a major soil moisture pulse occurs after spring thaw.
    %   See Supplementary Fig. 2 for details.
    %
    % The input format for start_datee, end_datee, linking_date, pf_datetime, final_thawed_date is as follows:
    %       start_datee / end_datee / linking_date / pf_datetime / final_thawed_date = 
    %                             [
    %                              datetime('2016-07-26');...
    %                              datetime('2017-07-24');...
    %                              datetime('2018-08-03');...
    %                              datetime('2019-08-02');...
    %                              datetime('2020-07-22');...
    %                              datetime('2021-07-07');...
    %                             ]; 
    % Each row corresponds to a specific year and defines the fitting start date, end date, or the sliding-linking date
    % used in the decomposition analysis.
    %
    % The input format for theta_pf is as follows:
    % theta_pf = [0.211;  0.208;  0.218;  0.199;  0.209;  0.213];
    % Each value corresponds to the pre-freezing soil moisture for a specific thaw-freeze cycle year.
    %
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i = 1:(len - 1)
        % Extract soil moisture data within the given time period
        idx_start = find(time_line == start_datee(i), 1);
        idx_end = find(time_line == end_datee(i), 1);
        if isempty(idx_start) || isempty(idx_end) || idx_start >= idx_end
            warning('Invalid start or end date in time_line, skipping year %d', years(i+1));
            continue;
        end
        raw_data = SM_check(idx_start:idx_end);
        x_vals = (1:length(raw_data))';
        x0 = x_vals(1);
        y0 = raw_data(1);

        % Fit with natural logarithmic function using fittype
        ft = fittype(@(a, b, x) a*log(b*x) + (y0 - a*log(b*x0)), ...
            'independent', 'x', 'coefficients', {'a','b'});
        log_fit = fit(x_vals, raw_data, ft, ...
            'StartPoint', [1, 1], 'Lower', [-Inf, 1e-5], 'Upper', [Inf, Inf]);

        y_fit = log_fit.a * log(log_fit.b * x_vals) + (y0 - log_fit.a * log(log_fit.b * x0));

        idx_link = find(time_line == linking_date(i), 1);
        if isempty(idx_link)
            warning('linking_date %s not found in time_line, skipping year %d', datestr(linking_date(i)), years(i+1));
            continue;
        end
        delta_filt = y_fit(1) - SM_check(idx_link);
        time_ss = linking_date(i) + caldays(0:length(y_fit)-1);
        theta_connect = y_fit - delta_filt;

        idx_ft = find(time_line == final_thawed_date(i), 1);
        idx_link_end = idx_link;
        if isempty(idx_ft)
            warning('final_thawed_date %s not found in time_line, skipping year %d', datestr(final_thawed_date(i)), years(i+1));
            continue;
        end
        theta1 = [SM_check(idx_ft:idx_link_end); theta_connect];
        figure
        plot(theta1)

        pf_end_date = datetime([num2str(years(i+1)) '-12-31']);
        claculate_time_line = pf_datetime(i):caldays(1):pf_end_date;
        theta_cut = SM_check(time_line >= pf_datetime(i) & time_line <= pf_end_date);

        theta_base = theta_crisis - theta_pf(i);
        theta_reward = theta1 - theta_crisis;
        theta_reward(theta_reward < 0) = 0;
        a = theta_reward ./ (theta_base + theta_reward);
        b = theta_base ./ (theta_base + theta_reward);

        time_index = final_thawed_date(i):caldays(1):time_ss(end);
        Started_base_time = datenum(final_thawed_date(i));
        idx_t1 = find(theta1 >= theta_crisis, 1, 'last');
        if isempty(idx_t1)
            warning('No point in theta1 exceeds threshold, skipping year %d', years(i+1));
            continue;
        end
        t1 = min([idx_t1, length(a), length(b), length(theta1)]);
        t1_date = time_index(t1);

        x_reward = (1:t1)';
        theta_reward_calcu = theta1(1:t1);
        E_reward = zeros(t1,1);
        for ii = 2:t1
            E_reward(ii) = trapz(x_reward(1:ii), max(0, theta_reward_calcu(1:ii) - theta_crisis));
        end

        theta_base_value = max(0, theta_crisis - theta_pf(i));
        E_base = theta_base_value * (x_reward - x_reward(1));
        E_base(1) = 0;
        E_SMR = a(1:t1) .* E_reward(1:t1) + b(1:t1) .* E_base(1:t1);

        above_threshold = SM_check > theta_crisis;
        excess_moisture = zeros(size(SM_check));
        excess_moisture(above_threshold) = SM_check(above_threshold) - theta_crisis;

        startt = datetime([num2str(years(i+1)) '-01-01']);
        endd = t1_date;
        time_mask = time_line >= startt & time_line <= endd;
        SM_calcu = excess_moisture(time_mask);
        times = time_line(time_mask);
        times_in_days = days(times - startt) + 1;
        E_all = trapz(times_in_days, SM_calcu);
        E_prec = E_all - E_reward(end);

        Results(i).year = years(i+1);
        Results(i).t1_date = t1_date;
        Results(i).theta_cut = theta_cut;
        Results(i).theta_connect = theta_connect;
        Results(i).theta1 = theta1;
        Results(i).theta_reward = theta_reward;
        Results(i).theta_base = theta_base;
        Results(i).E_reward = E_reward;
        Results(i).E_base = E_base;
        Results(i).E_SMR = E_SMR;
        Results(i).E_SMR_final = E_SMR(end);
        Results(i).E_all = E_all;
        Results(i).E_prec = E_prec;
        Results(i).t1 = t1;
        Results(i).time_index = time_index;
        Results(i).claculate_time_line = claculate_time_line;
    end
end
