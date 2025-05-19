% Synaptic parameters
g_STN_GPe = 0.1;  % Max conductance from STN to GPe (mS/cm^2)
E_STN_GPe = 0;    % Reversal potential for excitatory synapse (mV)

g_GPe_STN = 0.1;  % Max conductance from GPe to STN (mS/cm^2)
E_GPe_STN = -80;  % Reversal potential for inhibitory synapse (mV)

tau_syn = 5;      % Synaptic time constant (ms)

% Time parameters
dt = 0.01;                  % Time step (ms)
t = 0:dt:1000;              % Time vector (ms)
nSteps = length(t);


dw = 0.1; 
W_GPe2STN = 0.3:dw:2.5;
nWeight = length(W_GPe2STN);

% Initialize membrane potentials
V_STN = zeros(1, nSteps);
V_GPe = zeros(1, nSteps);
V_STN(1) = -65;  % Initial membrane potential for STN (mV)
V_GPe(1) = -65;  % Initial membrane potential for GPe (mV)

% Initialize synaptic gating variables
s_STN_GPe = zeros(1, nSteps);
s_GPe_STN = zeros(1, nSteps);

% Initialize calcium concentrations
Ca_STN = zeros(1, nSteps);
Ca_GPe = zeros(1, nSteps);
Ca_STN(1) = 1e-4;
Ca_GPe(1) = 1e-4;

% Initialize injected currents
Iinj_STN = zeros(1, nSteps);
Iinj_GPe = zeros(1, nSteps);

function I = generate_beta_poisson_bursts(t, beta_freq, base_rate, burst_rate, amp, jitter)
    dt = t(2) - t(1);
    nSteps = length(t);
    I = zeros(1, nSteps);

    osc = sin(2*pi*beta_freq*t/1000);  % beta envelope
    modulated_rate = base_rate + (burst_rate - base_rate) * (osc > 0);  % half-wave rectified

    for i = 1:nSteps
        if rand < modulated_rate(i) * dt / 1000
            pulse_width_steps = round(jitter / dt);
            end_idx = min(i + pulse_width_steps - 1, nSteps);
            I(i:end_idx) = I(i:end_idx) + amp;
        end
    end
end

beta_freq = 20;              % Hz
base_rate = 5;               % baseline (Hz)
burst_rate = 300;            % in-burst rate (Hz)
amp = 55;                    % current pulse amplitude (μA/cm²)
jitter = 0.5;                  % ms per spike

Ictx = generate_beta_poisson_bursts(t, beta_freq, base_rate, burst_rate, amp, jitter);
V_STN_all = zeros(nWeight, nSteps);
V_GPe_all = zeros(nWeight, nSteps);
phase_diff_deg = zeros(1, nWeight);

for j = 1:nWeight
    % Reinitialize for each weight
    V_STN = zeros(1, nSteps); V_STN(1) = -65;
    V_GPe = zeros(1, nSteps); V_GPe(1) = -65;
    Ca_STN = zeros(1, nSteps); Ca_STN(1) = 1e-4;
    Ca_GPe = zeros(1, nSteps); Ca_GPe(1) = 1e-4;
    s_STN_GPe = zeros(1, nSteps);
    s_GPe_STN = zeros(1, nSteps);
    Iinj_STN = zeros(1, nSteps);
    Iinj_GPe = zeros(1, nSteps);

    state_STN = struct('m', 0, 'h', 1, 'n', 0, 'r', 0, 's', 1, 'f', 0);
    state_GPe = struct( ...
        'm', 1 / (1 + exp(-( -65 + 35)/7)), ...
        'h', 1 / (1 + exp(( -65 + 60)/7)), ...
        'n', 1 / (1 + exp(-( -65 + 30)/10)), ...
        'q', 1 / (1 + exp(-( -65 + 25)/5)), ...
        'f', 1 / (1 + exp(( -65 + 75)/5.5)) ...
    );

    for i = 2:nSteps
        I_syn_STN = g_GPe_STN * s_GPe_STN(i-1) * (V_STN(i-1) - E_GPe_STN);
        I_syn_GPe = g_STN_GPe * s_STN_GPe(i-1) * (V_GPe(i-1) - E_STN_GPe);
        
        Iinj_STN(i-1) = Ictx(i-1) - W_GPe2STN(j) * V_GPe(i-1) + I_syn_STN;
        Iinj_GPe(i-1) = I_syn_GPe;

        [V_STN(i), ~, Ca_STN(i), state_STN] = STN_Reduced_Step(V_STN(i-1), Ca_STN(i-1), Iinj_STN(i-1), dt, state_STN);
        [V_GPe(i), ~, Ca_GPe(i), state_GPe] = GPe_Reduced_Step(V_GPe(i-1), Ca_GPe(i-1), Iinj_GPe(i-1), dt, state_GPe);

        % Synaptic gating update
        s_STN_GPe(i) = s_STN_GPe(i-1) + dt * (-s_STN_GPe(i-1)/tau_syn);
        s_GPe_STN(i) = s_GPe_STN(i-1) + dt * (-s_GPe_STN(i-1)/tau_syn);

        if V_STN(i-1) > 0 && V_STN(i-2) <= 0
            s_STN_GPe(i) = s_STN_GPe(i) + 1;
        end
        if V_GPe(i-1) > 0 && V_GPe(i-2) <= 0
            s_GPe_STN(i) = s_GPe_STN(i) + 1;
        end
    end

    V_STN_all(j, :) = V_STN;
    V_GPe_all(j, :) = V_GPe;
end


t = 0:0.01:1000;  % ms

phase_diff_STN_GPe_range(t, V_STN_all, V_GPe_all, W_GPe2STN, beta_freq);

figure;
plot(t, V_STN, 'b', t, V_GPe, 'r');
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
legend('STN', 'GPe');
title('STN-GPe Network Dynamics');