function [V_new, I_total, Ca_new, state] = STN_Reduced_Step(V, Ca, Iinj, dt, state)
    % Reversal potentials (mV)
    ENa = 55; EK = -90; ECa = 120; EHCN = -40; EL = -65;

    % Max conductances (mS/cm²)
    gNaF = 60; gKDR = 20; gCaT = 1.8;
    gSK = 2.0; gHCN = 0.4; gL = 0.1;

    % Unpack gating variables
    m = state.m; h = state.h; n = state.n;
    r = state.r; s = state.s; f = state.f;

    % NaF gating variables
    a_m = 0.1*(V+40)/(1 - expc(-(V+40)/10));
    b_m = 4 * expc(-(V+65)/18);
    m_inf = a_m / (a_m + b_m); tau_m = 1 / (a_m + b_m);

    a_h = 0.07 * expc(-(V+65)/20);
    b_h = 1 / (1 + expc(-(V+35)/10));
    h_inf = a_h / (a_h + b_h); tau_h = 1 / (a_h + b_h);

    % KDR gating variables
    a_n = 0.01*(V+55)/(1 - expc(-(V+55)/10));
    b_n = 0.125*expc(-(V+65)/80);
    n_inf = a_n / (a_n + b_n); tau_n = 1 / (a_n + b_n);

    % CaT gating variables
    r_inf = 1 / (1 + expc(-(V+60)/7.5)); tau_r = 5;
    s_inf = 1 / (1 + expc((V+70)/5));  tau_s = 30;

    % HCN gating variables
    f_inf = 1 / (1 + expc((V + 75)/5.5)); tau_f = 300;

    % SK gating variable
    w = 1 / (1 + (0.3 / max(Ca, 1e-9))^4);

    % Ionic currents
    INa = gNaF * m^3 * h * (V - ENa);
    IK = gKDR * n^4 * (V - EK);
    ICaT = gCaT * r^3 * s * (V - ECa);
    IsK = gSK * w * (V - EK);
    IHCN = gHCN * f * (V - EHCN);
    IL = gL * (V - EL);

    I_total = INa + IK + ICaT + IsK + IHCN + IL;
    dV = Iinj - I_total;
    V_new = V + dt * dV;

    % Calcium dynamics
    dCa = -0.01 * ICaT - (Ca - 1e-4)/80;
    Ca_new = max(Ca + dt * dCa, 1e-7);

    % Update gating variables
    state.m = clamp(m + dt * (m_inf - m)/tau_m);
    state.h = clamp(h + dt * (h_inf - h)/tau_h);
    state.n = clamp(n + dt * (n_inf - n)/tau_n);
    state.r = clamp(r + dt * (r_inf - r)/tau_r);
    state.s = clamp(s + dt * (s_inf - s)/tau_s);
    state.f = clamp(f + dt * (f_inf - f)/tau_f);
end


function y = clamp(x)
    y = min(max(x, 0), 1);
end

function y = expc(x)
    y = exp(min(max(x, -50), 50));
end