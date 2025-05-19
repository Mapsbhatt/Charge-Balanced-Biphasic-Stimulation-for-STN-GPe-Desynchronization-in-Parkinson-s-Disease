function [V_new, I_total, Ca_new, state] = GPe_Reduced_Step(V, Ca, Iinj, dt, state)
    % Reversal potentials (mV)
    ENa = 55; EK = -90; ECa = 120; EHCN = -40; EL = -65;

    % Max conductances (mS/cm²)
    gNa = 60; gK = 25;
    gCa = 0.5; gSK = 1.5;
    gHCN = 0.2; gL = 0.1;
    gNaP = 0.4;  % Persistent Na+
    

    % Unpack gating variables
    m = state.m; h = state.h; n = state.n;
    q = state.q; f = state.f;

    % INaP gating (steady-state)
    mNaP_inf = 1 / (1 + exp(-(V + 52)/5.3));
    INaP = gNaP * mNaP_inf * (V - ENa);

    % Na gating variables
    m_inf = 1 / (1 + exp(-(V+35)/7));
    h_inf = 1 / (1 + exp((V+60)/7));
    tau_m = 0.3; tau_h = 1.0;

    % K gating variables
    n_inf = 1 / (1 + exp(-(V+30)/10));
    tau_n = 1.0;

    % Ca gating variables
    q_inf = 1 / (1 + exp(-(V+25)/5));
    tau_q = 5;

    % HCN gating variables
    f_inf = 1 / (1 + exp((V + 75)/5.5));
    tau_f = 300;

    % SK gating variable
    w = 1 / (1 + (0.3 / max(Ca, 1e-9))^4);

    % Ionic currents
    INa = gNa * m^3 * h * (V - ENa);
    IK = gK * n^4 * (V - EK);
    ICa = gCa * q^2 * (V - ECa);
    IsK = gSK * w * (V - EK);
    IHCN = gHCN * f * (V - EHCN);
    IL = gL * (V - EL);

    I_total = INa + INaP + IK + ICa + IsK + IHCN + IL;
    dV = Iinj - I_total;
    V_new = V + dt * dV;

    % Calcium dynamics
    dCa = -0.01 * ICa - (Ca - 1e-4)/80;
    Ca_new = max(Ca + dt * dCa, 1e-7);

    % Update gating variables
    state.m = clamp(m + dt * (m_inf - m)/tau_m);
    state.h = clamp(h + dt * (h_inf - h)/tau_h);
    state.n = clamp(n + dt * (n_inf - n)/tau_n);
    state.q = clamp(q + dt * (q_inf - q)/tau_q);
    state.f = clamp(f + dt * (f_inf - f)/tau_f);
end


function y = clamp(x)
    y = min(max(x, 0), 1);
end

function y = expc(x)
    y = exp(min(max(x, -50), 50));
end