function x_next = olivo_model(x, u)
    % =====================================================================
    % En NLMPC, el vector 'u' agrupa TODAS las entradas (MV y MD)
    % u(1): Riego I(k) [MV]
    % u(2): Precipitacion [MD]
    % u(3): ET0 [MD]
    % u(4): Velocidad del viento u2 [MD]
    % u(5): Humedad relativa mínima hrmin [MD]
    % u(6): Kcb_k [MD]
    % =====================================================================

    % 1. Desempaquetar variables 
    R_k     = x(1);
    I_k     = u(1);
    Precip  = u(2);
    et0_k   = u(3);
    u2_k    = u(4);
    hrmin_k = u(5);
    Kcb_k   = u(6);
    
    % 2. Parámetros fijos de la Planta
    Rfc = 168;
    h = 3; %antes estaba en 4 m porque 3 - 5m pero al ser olivar intensivo mejor 3m
    kcmin = 0.15;
    do = 1;
    a = 0.03;
    umbral_lluvia = 3.5;
    fwgoteo = 0.3;

    % 3. Aplicamos Ke y Kcmax
    Kcmax_k = max(1.2 + (0.04*(u2_k - 2) - 0.004*(hrmin_k - 45)) * (h/3)^0.3, Kcb_k + 0.05);
    fc_k = ((Kcb_k - kcmin) / (Kcmax_k - kcmin)) ^ (1 + 0.5*h);
    
    % 4. Saturación de fc
    if fc_k < 0
        fc_k = 0;
    elseif fc_k > 0.99
        fc_k = 0.99;
    end

    % 5. Dinámica de suelo mojado
    if Precip > umbral_lluvia 
        fw_k = 1.0; 
    elseif I_k > 0
        fw_k = fwgoteo*(1-((2/3)*fc_k)); 
    else 
        fw_k = fwgoteo*(1-((2/3)*fc_k));
    end

    few_k = min(1 - fc_k, fw_k);
    Ke_k = few_k * Kcmax_k;

    % 6. Evapotranspiración y Drenaje
    ET_k = (Kcb_k + Ke_k) * et0_k;
    DP_k = do * ((exp(a * R_k) - 1) / (exp(a * Rfc) - 1));

    % 7. Ecuación final del Balance de Masas
    x_next = R_k + I_k + Precip - ET_k - DP_k;
end