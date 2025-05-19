% Revised STN-GPe simulation with heterogeneous external inputs

% Synaptic parameters
g_STN_GPe = 0.1;  % STN to GPe excitatory
E_STN_GPe = 0;    % mV

g_GPe_STN = 0.1;  % GPe to STN inhibitory
E_GPe_STN = -80;  % mV

tau_syn = 5;      % ms

% Time parameters
dt = 0.01;
t = 0:dt:1000;  % ms
nSteps = length(t);

% Membrane potentials
V_STN = zeros(1, nSteps);
V_GPe = zeros(1, nSteps);
V_STN(1) = -65;
V_GPe(1) = -65;

% Synaptic gating
s_STN_GPe = zeros(1, nSteps);
s_GPe_STN = zeros(1, nSteps);

% Calcium
Ca_STN = zeros(1, nSteps);
Ca_GPe = zeros(1, nSteps);
Ca_STN(1) = 1e-4;
Ca_GPe(1) = 1e-4;

% Injected currents
Iinj_STN = zeros(1, nSteps);
Iinj_GPe = zeros(1, nSteps);

% Initial states
state_STN = struct('m',0,'h',1,'n',0,'r',0,'s',1,'f',0);
state_GPe = struct('m',1/(1+exp(-( -65 + 35)/7)), 'h',1/(1+exp(( -65 + 60)/7)), 'n',1/(1+exp(-( -65 + 30)/10)), 'q',1/(1+exp(-( -65 + 25)/5)), 'f',1/(1+exp(( -65 + 75)/5.5)));

% Heterogeneous beta Poisson cortical inputs
beta_freq = 20;
base_rate = 5;
burst_rate = 300;
amp = 55;

%% Change the jitter 
% jitter = 0.5; 
jitter = 0.8;

gap = 0.5;  % ms inter-phase gap
Ictx_STN = generate_beta_poisson_biphasic(t, beta_freq, base_rate, burst_rate, amp, jitter, gap);
Ictx_GPe = generate_beta_poisson_biphasic(t + 25, beta_freq, base_rate, burst_rate, amp, jitter, gap); %25ms delay


% Simulation
for i = 2:nSteps-1
    I_syn_STN = g_GPe_STN * s_GPe_STN(i) * (V_STN(i) - E_GPe_STN);
    I_syn_GPe = g_STN_GPe * s_STN_GPe(i) * (V_GPe(i) - E_STN_GPe);
    
    Iinj_STN(i) = Ictx_STN(i) + I_syn_STN;   % Stimulation 
    Iinj_GPe(i) = Ictx_GPe(i) + I_syn_GPe;

    [V_STN(i+1), ~, Ca_STN(i+1), state_STN] = STN_Reduced_Step(V_STN(i), Ca_STN(i), Iinj_STN(i), dt, state_STN);
    [V_GPe(i+1), ~, Ca_GPe(i+1), state_GPe] = GPe_Reduced_Step(V_GPe(i), Ca_GPe(i), Iinj_GPe(i), dt, state_GPe);

    s_STN_GPe(i+1) = s_STN_GPe(i) + dt * (-s_STN_GPe(i)/tau_syn);
    s_GPe_STN(i+1) = s_GPe_STN(i) + dt * (-s_GPe_STN(i)/tau_syn);

    if V_STN(i) > 0 && V_STN(i-1) <= 0
        s_STN_GPe(i+1) = s_STN_GPe(i+1) + 1;
    end
    if V_GPe(i) > 0 && V_GPe(i-1) <= 0
        s_GPe_STN(i+1) = s_GPe_STN(i+1) + 1;
    end
end

% Plot
figure;
plot(t, V_STN, 'b', t, V_GPe, 'r');
xlabel('Time (ms)'); ylabel('Membrane Potential (mV)');
legend('STN', 'GPe');
title('STN-GPe Dynamics with Heterogeneous Cortical Inputs');

% Supporting functions

function I = generate_beta_poisson_biphasic(t, beta_freq, base_rate, burst_rate, amp, jitter, gap)
    dt = t(2) - t(1);
    nSteps = length(t);
    I = zeros(1, nSteps);

    osc = sin(2*pi*beta_freq*t/1000);  % beta envelope
    modulated_rate = base_rate + (burst_rate - base_rate) * (osc > 0);  % rectified beta modulation

    pulse_width_steps = round(jitter / dt);       % duration of each phase
    gap_steps = round(gap / dt);                  % inter-phase gap
    total_steps = 2 * pulse_width_steps + gap_steps;

    for i = 1:nSteps
        if rand < modulated_rate(i) * dt / 1000
            idx_end = min(i + total_steps - 1, nSteps);
            % Build biphasic pulse: [-A, gap, +A]
            biphasic = [-amp * ones(1, pulse_width_steps), ...
                         zeros(1, gap_steps), ...
                         amp * ones(1, pulse_width_steps)];
            len = min(length(biphasic), idx_end - i + 1);
            I(i:i+len-1) = I(i:i+len-1) + biphasic(1:len);
        end
    end


end
