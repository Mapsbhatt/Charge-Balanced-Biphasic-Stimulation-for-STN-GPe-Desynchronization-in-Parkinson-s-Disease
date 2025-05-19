function phase_diff_rad = phase_diff_STN_GPe_spikebased(t, V_STN_all, V_GPe_all, W_values, beta_freq)
% Calculates the spike-based phase difference between STN and GPe neurons
% within 550–650 ms, based on nearest-matching spike peaks

nW = length(W_values);
phase_diff_rad = NaN(1, nW);
dt = t(2) - t(1);
T_beta = 1000 / beta_freq;  % beta period in ms

for i = 1:nW
    v_stn = V_STN_all(i, :);
    v_gpe = V_GPe_all(i, :);

    % Window: 550–650 ms
    idx_range = find(t >= 550 & t <= 650);
    t_win = t(idx_range);
    v_stn_win = v_stn(idx_range);
    v_gpe_win = v_gpe(idx_range);

    % Find peaks above 5 mV
    [pks_stn, locs_stn] = findpeaks(v_stn_win, 'MinPeakHeight', 5);
    [pks_gpe, locs_gpe] = findpeaks(v_gpe_win, 'MinPeakHeight', 5);

    if isempty(locs_stn) || isempty(locs_gpe)
        warning('No peaks found for W = %.2f', W_values(i));
        continue;
    end

    % Use the largest STN peak (assume it's most reliable)
    [~, idx_stn_max] = max(pks_stn);
    t_peak_stn = t_win(locs_stn(idx_stn_max));

    % Find closest GPe peak
    [~, idx_gpe_closest] = min(abs(t_win(locs_gpe) - t_peak_stn));
    t_peak_gpe = t_win(locs_gpe(idx_gpe_closest));

    % Calculate phase difference as fraction of beta cycle
    delta_t = t_peak_stn - t_peak_gpe;
    phi = 2 * pi * delta_t / T_beta;
    phase_diff_rad(i) = wrapToPi(phi);
end

% Plot
figure;
plot(W_values, phase_diff_rad, '-o', 'LineWidth', 2);
xlabel('W_{GPe→STN}');
ylabel('Phase Difference (radians)');
title('STN-GPe Phase Lag (Spike-Based, 550–650 ms)');
yticks([-pi -pi/2 0 pi/2 pi]);
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
grid on;
end