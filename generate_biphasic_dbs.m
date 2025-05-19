function I_dbs = generate_biphasic_dbs(t, pulse_freq, pulse_amp, pulse_width)
    % Generates a charge-balanced biphasic square pulse train.
    %
    % Parameters:
    % t           - Time vector (ms)
    % pulse_freq  - Pulse frequency (Hz)
    % pulse_amp   - Pulse amplitude (μA/cm²)
    % pulse_width - Width of each phase (ms)
    %
    % Returns:
    % I_dbs       - Biphasic current waveform (μA/cm²)

    dt = t(2) - t(1);  % Time step (ms)
    I_dbs = zeros(size(t));
    period = 1000 / pulse_freq;  % Pulse period in ms
    phase_steps = round(pulse_width / dt);

    for i = 1:length(t)
        cycle_pos = mod(t(i), period);
        if cycle_pos < pulse_width
            I_dbs(i) = -pulse_amp;
        elseif cycle_pos < 2 * pulse_width
            I_dbs(i) = pulse_amp;
        else
            I_dbs(i) = 0;
        end
    end
end