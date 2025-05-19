function [V, Iion, Iinj, Ca] = GPe_Reduced(t, Iinj)
% Reduced GPe neuron model with INa, IK, ICa, IsKCa, IHCN, IL

dt = t(2) - t(1);
nSteps = length(t);
V = zeros(1, nSteps); V(1) = -65;
Ca = zeros(1, nSteps); Ca(1) = 1e-4;
Iion = zeros(1, nSteps);

% Reversal potentials (mV)
ENa = 55; EK = -90; ECa = 120; EHCN = -40; EL = -65;

% Max conductances (mS/cm²)
gNa = 60; gK = 25;
gCa = 0.5; gSK = 1.5;
gHCN = 0.2; gL = 0.1;
gNaP = 0.4;  % mS/cm², adjust to get tonic firing

% Initial gate values (steady-state at -65 mV)
V0 = V(1); Ca0 = Ca(1);
m = 1 / (1 + exp(-(V0+35)/7));
h = 1 / (1 + exp((V0+60)/7));
n = 1 / (1 + exp(-(V0+30)/10));
q = 1 / (1 + exp(-(V0+25)/5));  % Ca gate
f = 1 / (1 + exp((V0 + 75)/5.5));  % HCN

for i = 1:nSteps-1
    Vnow = V(i); Ca_now = Ca(i);
    % INaP gating (steady-state only)
    mNaP_inf = 1 / (1 + exp(-(Vnow + 52)/5.3));  % from Koelman & Lowery
    INaP = gNaP * mNaP_inf * (Vnow - ENa);

    % Na activation/inactivation
    m_inf = 1 / (1 + exp(-(Vnow+35)/7));
    h_inf = 1 / (1 + exp((Vnow+60)/7));
    tau_m = 0.3; tau_h = 1.0;
    dm = (m_inf - m)/tau_m;
    dh = (h_inf - h)/tau_h;

    % K activation
    n_inf = 1 / (1 + exp(-(Vnow+30)/10));
    tau_n = 1.0;
    dn = (n_inf - n)/tau_n;

    % HVA Ca²⁺
    q_inf = 1 / (1 + exp(-(Vnow+25)/5));
    tau_q = 5;
    dq = (q_inf - q)/tau_q;

    % HCN
    f_inf = 1 / (1 + exp((Vnow + 75)/5.5));
    tau_f = 300;
    df = (f_inf - f)/tau_f;

    % SK gating (Ca-dependent Hill function)
    w = 1 / (1 + (0.3 / max(Ca_now, 1e-9))^4);

    % Ionic currents
    INa = gNa * m^3 * h * (Vnow - ENa);
    IK = gK * n^4 * (Vnow - EK);
    ICa = gCa * q^2 * (Vnow - ECa);
    IsK = gSK * w * (Vnow - EK);
    IHCN = gHCN * f * (Vnow - EHCN);
    IL = gL * (Vnow - EL);

    I_total = INa + INaP + IK + ICa + IsK + IHCN + IL;
    dV = Iinj(i) - I_total;
    V(i+1) = Vnow + dt * dV;

    % Calcium dynamics
    dCa = -0.01 * ICa - (Ca_now - 1e-4)/80;
    Ca(i+1) = max(Ca_now + dt * dCa, 1e-7);

    % Update gates
    m = clamp(m + dt*dm); h = clamp(h + dt*dh);
    n = clamp(n + dt*dn);
    q = clamp(q + dt*dq);
    f = clamp(f + dt*df);

    Iion(i) = I_total;
end
end

function y = clamp(x)
y = min(max(x, 0), 1);
end


t = 0:0.01:1000;
Iinj = zeros(size(t));
% Iinj(20001:20020) = 0;  % Now it's a 0.2 ms pulse at 200 ms

[V, Iion, Iinj, Ca] = GPe_Reduced(t, Iinj);
plot(t, V);
xlabel('Time (ms)'); ylabel('Membrane Potential (mV)');
title('Reduced GPe Neuron Model');