function I = generate_beta_poisson_bursts(t, beta_freq, base_rate, burst_rate, amp, jitter)
% Generate Poisson spike bursts modulated by beta-band oscillations
%
% Inputs:
%   t           = time vector (ms)
%   beta_freq   = beta modulation frequency (Hz), e.g. 20
%   base_rate   = base firing rate (Hz) outside burst (e.g. 5)
%   burst_rate  = firing rate inside burst (e.g. 300 Hz)
%   amp         = amplitude of each spike (μA/cm²)
%   jitter      = pulse width (ms) for each spike
%
% Output:
%   I = current injection array (same size as t)

dt = t(2) - t(1);
nSteps = length(t);
I = zeros(1, nSteps);

% Create a sine wave envelope for beta modulation
osc = sin(2*pi*beta_freq*t/1000);  % beta envelope (unit amplitude)
modulated_rate = base_rate + (burst_rate - base_rate) * (osc > 0);  % half-wave rectified

% Generate Poisson spikes with modulated rate
for i = 1:nSteps
    if rand < modulated_rate(i) * dt / 1000  % firing probability at this time
        pulse_width_steps = round(jitter / dt);
        end_idx = min(i + pulse_width_steps - 1, nSteps);
        I(i:end_idx) = I(i:end_idx) + amp;
    end
end

end


