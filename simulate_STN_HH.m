function [V, Iion, Iinj, Ca] = STN_Reduced(t, Iinj)
% Reduced STN neuron model with spike, burst, HCN dynamics

dt = t(2) - t(1); nSteps = length(t);
V = zeros(1, nSteps); V(1) = -65;
Ca = zeros(1, nSteps); Ca(1) = 1e-4;
%Iinj = zeros(1, nSteps); Iion = zeros(1, nSteps);

% Reversal potentials
ENa = 55; EK = -90; ECa = 120; EHCN = -40; EL = -65;

% Conductances (mS/cm²)
gNaF = 60; gKDR = 20; gCaT = 1.8;
gSK = 2.0; gHCN = 0.4; gL = 0.1;

% === Initialize gating variables at steady-state ===
V0 = V(1); Ca0 = Ca(1);

% NaF
a_m = 0.1*(V0+40)/(1 - expc(-(V0+40)/10));
b_m = 4 * expc(-(V0+65)/18);
m = a_m / (a_m + b_m);

a_h = 0.07 * expc(-(V0+65)/20);
b_h = 1 / (1 + expc(-(V0+35)/10));
h = a_h / (a_h + b_h);

% KDR
a_n = 0.01*(V0+55)/(1 - expc(-(V0+55)/10));
b_n = 0.125*expc(-(V0+65)/80);
n = a_n / (a_n + b_n);

% CaT
r = 1 / (1 + expc(-(V0+60)/7.5));
s = 1 / (1 + expc((V0+70)/5));

% HCN
f = 1 / (1 + expc((V0 + 75)/5.5));

for i = 1:nSteps-1
    Vnow = V(i); Ca_now = Ca(i);

    % Pulse input
    %Iinj(i) = (t(i) >= pd && t(i) <= pd + pw) * pa;
    % fc = 10;  % Hz
    % Iinj = 50 * (sin(2*pi*fc*t/1000) > 0.96);  % 10 Hz sine wave envelope
    
  
    % NaF gating
    a_m = 0.1*(Vnow+40)/(1 - expc(-(Vnow+40)/10));
    b_m = 4 * expc(-(Vnow+65)/18);
    m_inf = a_m / (a_m + b_m); tau_m = 1 / (a_m + b_m);

    a_h = 0.07 * expc(-(Vnow+65)/20);
    b_h = 1 / (1 + expc(-(Vnow+35)/10));
    h_inf = a_h / (a_h + b_h); tau_h = 1 / (a_h + b_h);

    dm = (m_inf - m)/tau_m;
    dh = (h_inf - h)/tau_h;

    % KDR
    a_n = 0.01*(Vnow+55)/(1 - expc(-(Vnow+55)/10));
    b_n = 0.125*expc(-(Vnow+65)/80);
    n_inf = a_n / (a_n + b_n); tau_n = 1 / (a_n + b_n);
    dn = (n_inf - n)/tau_n;

    % CaT
    r_inf = 1 / (1 + expc(-(Vnow+60)/7.5)); tau_r = 5;
    s_inf = 1 / (1 + expc((Vnow+70)/5)); tau_s = 30;
    dr = (r_inf - r)/tau_r; ds = (s_inf - s)/tau_s;

    % HCN
    f_inf = 1 / (1 + expc((Vnow + 75)/5.5)); tau_f = 300;
    df = (f_inf - f)/tau_f;

    % SK gating (Hill function)
    w = 1 / (1 + (0.3 / max(Ca_now, 1e-9))^4);

    % Currents
    INa = gNaF * m^3 * h * (Vnow - ENa);
    IK = gKDR * n^4 * (Vnow - EK);
    ICaT = gCaT * r^3 * s * (Vnow - ECa);
    IsK = gSK * w * (Vnow - EK);
    IHCN = gHCN * f * (Vnow - EHCN);
    IL = gL * (Vnow - EL);

    I_total = INa + IK + ICaT + IsK + IHCN + IL;
    dV = Iinj(i) - I_total;
    V(i+1) = Vnow + dt * dV;

    % Calcium dynamics
    dCa = -0.01 * ICaT - (Ca_now - 1e-4)/80;
    Ca(i+1) = max(Ca_now + dt * dCa, 1e-7);

    % Update gates
    m = clamp(m + dt * dm); h = clamp(h + dt * dh);
    n = clamp(n + dt * dn);
    r = clamp(r + dt * dr); s = clamp(s + dt * ds);
    f = clamp(f + dt * df);

    Iion(i) = I_total;
end
end

function y = expc(x)
y = exp(min(max(x, -50), 50));
end

function y = clamp(x)
y = min(max(x, 0), 1);
end

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


t = 0:0.01:1000;  % ms
beta_freq = 20;              % Hz
base_rate = 5;               % baseline (Hz)
burst_rate = 300;            % in-burst rate (Hz)
amp = 55;                    % current pulse amplitude (μA/cm²)
jitter = 0.5;                  % ms per spike

Iinj = generate_beta_poisson_bursts(t, beta_freq, base_rate, burst_rate, amp, jitter);

[V, Iion, Iinj, Ca] = STN_Reduced(t, Iinj);
figure;
plot(t, V); xlabel('Time (ms)'); ylabel('Membrane Potential (mV)');
title('Reduced STN Neuron: Spike + Burst + HCN Dynamics');


