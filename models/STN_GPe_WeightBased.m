% Synaptic parameters
g_STN_GPe = 0.1;  % Max conductance from STN to GPe (mS/cm^2)
E_STN_GPe = 0;    % Reversal potential for excitatory synapse (mV)

  % Max conductance from GPe to STN (mS/cm^2)
E_GPe_STN = -80;  % Reversal potential for inhibitory synapse (mV)

tau_syn = 5;      % Synaptic time constant (ms)

% Time parameters
dt = 0.01;                  % Time step (ms)
t = 0:dt:1000;              % Time vector (ms)
nSteps = length(t);


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
I_dbs = zeros(1, nSteps);

state_STN = struct( ...
    'm', 0, ...
    'h', 1, ...
    'n', 0, ...
    'r', 0, ...
    's', 1, ...
    'f', 0 ...
);

state_GPe = struct( ...
    'm', 1 / (1 + exp(-( -65 + 35)/7)), ...
    'h', 1 / (1 + exp(( -65 + 60)/7)), ...
    'n', 1 / (1 + exp(-( -65 + 30)/10)), ...
    'q', 1 / (1 + exp(-( -65 + 25)/5)), ...
    'f', 1 / (1 + exp(( -65 + 75)/5.5)) ...
);

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


% DBS pulse parameters

scale = 0.9;         % Matching your W_GPe2STN offset
pw = 3;              % Each phase duration in ms

% Generate biphasic pulse train
g_GPe_STN = 0.1; 

%% Vary the value from 0.4 to 1.5 to see in sync action
W_GPe2STN = 0.4; % Normal functioning
%W_GPe2STN = 1.5; % Parkinson Condition

for i = 1:nSteps-1
        % Compute synaptic currents
        I_syn_STN = g_GPe_STN * s_GPe_STN(i) * (V_STN(i) - E_GPe_STN);
        I_syn_GPe = g_STN_GPe * s_STN_GPe(i) * (V_GPe(i) - E_STN_GPe);
        I_dbs(i) = scale * V_GPe(i) * ((mod(floor(t(i)/(2*pw)), 2) == 0) - 0.5) * 2;
        
        % Total injected currents including synaptic inputs
        Iinj_STN(i) = I_syn_STN;
        Iinj_STN(i) = Ictx(i) - W_GPe2STN * V_GPe(i); 
        %Iinj_STN(i) = Ictx(i) - W_GPe2STN * V_GPe(i) + I_dbs(i);  % Excitatory cortex + inhibitory GPe
        Iinj_GPe(i) = I_syn_GPe;
        
        % Update STN neuron
        [V_STN(i+1), ~, Ca_STN(i+1), state_STN] = STN_Reduced_Step(V_STN(i), Ca_STN(i), Iinj_STN(i), dt, state_STN);
        
        % Update GPe neuron
        [V_GPe(i+1), ~, Ca_GPe(i+1), state_GPe] = GPe_Reduced_Step(V_GPe(i), Ca_GPe(i), Iinj_GPe(i), dt, state_GPe);
        
        % Update synaptic gating variables
        s_STN_GPe(i+1) = s_STN_GPe(i) + dt * (-s_STN_GPe(i)/tau_syn);
        s_GPe_STN(i+1) = s_GPe_STN(i) + dt * (-s_GPe_STN(i)/tau_syn);
        
        % Spike detection and synaptic gating update
        if V_STN(i) > 0 && V_STN(i-1) <= 0
            s_STN_GPe(i+1) = s_STN_GPe(i+1) + 1;
        end
        if V_GPe(i) > 0 && V_GPe(i-1) <= 0
            s_GPe_STN(i+1) = s_GPe_STN(i+1) + 1;
        end
end

figure;
% Plot STN Membrane Potential
yyaxis left;
plot(t, V_STN, 'b', 'LineWidth', 1.2);
ylabel('STN Membrane Potential (mV)');
ylim([-90 40]);

% Plot DBS Current
yyaxis right;
plot(t, I_dbs, 'r', 'LineWidth', 1.2);
ylabel('Injected DBS Current (μA/cm²)');
ylim([-50 50]);

xlabel('Time (ms)');
title('STN Voltage and DBS Current');
legend('STN V_m', 'DBS Current');
grid on;

t = 0:0.01:1000;  % ms



figure;
plot(t, V_STN, 'b', t, V_GPe, 'r');
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
legend('STN', 'GPe');
title('STN-GPe Network Dynamics');