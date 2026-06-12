% =========================================================================
% SCRIPT AISLADO: RUIDO DE PROCESO (OE5) - VERSIÓN BLINDADA
% =========================================================================
clear, clc; close all;

niveles_proceso = [1, 2, 3, 4, 5]; % Multiplicador del ruido de proceso (mm/día)
N_simulaciones = 30; 

Agua_Media = zeros(1, length(niveles_proceso));
Estres_Medio = zeros(1, length(niveles_proceso));

datos = readtable('Osuna3.csv'); 
Precipitacion = flipud(datos.Se11Precip);
et0 = flipud(datos.Se11ETo); 
u2 = flipud(datos.Se11VelViento); 
hrmin = flipud(datos.Se11HumMin);

N_dias = 365; theta_fc = 0.34; theta_wp = 0.20; theta_sat = 0.45; Zr = 1.20; h=3; kcmin=0.15; 
caudal_mm = 16/12; num_max_horas_riego = 8; 
Rfc = 1000 * (theta_fc - theta_wp) * Zr; 
Rmax = 1000 * (theta_sat - theta_wp) * Zr; 
R_min = Rfc * (1 - 0.65); 
R_objetivo = (Rfc + R_min) / 2; 

do = 1; a = 0.03; Qm = 0.1; Rm = 1; Np = 7; Nc = 7; 
Lluvia_extendida = [Precipitacion; Precipitacion(end-Np+1 : end)];
ET0_extendida = [et0; et0(end-Np+1 : end)]; 
u2_extendida = [u2; u2(end-Np+1 : end)];
hrmin_extendida = [hrmin; hrmin(end-Np+1 : end)];

Kcb = zeros(1, N_dias);    
for k = 1:N_dias
    if k < 60 || k >= 330
        Kcb_tab = 0.4;
    elseif k >= 60 && k < 90
        Kcb_tab = 0.55; 
    elseif k >= 90 && k < 180
        Kcb_tab = 0.55 + ((0.65 - 0.55) / 90) * (k - 90); 
    else
        Kcb_tab = 0.65; 
    end
    Kcb(k) = Kcb_tab + (0.04*(u2(k)-2) - 0.004*(hrmin(k)-45)) * (h/3)^0.3;
end
Kcb_extendida = [Kcb, 0.4 * ones(1, Np)];

nx = 1; ny = 1; 
nlobj = nlmpc(nx, ny, 'MV', 1, 'MD', [2 3 4 5 6]);
nlobj.Ts = 1; 
nlobj.PredictionHorizon = Np; 
nlobj.ControlHorizon = Nc;    
nlobj.Model.StateFcn = "olivo_model"; 
nlobj.Model.IsContinuousTime = false; 

nlobj.ManipulatedVariables(1).Min = 0; 
nlobj.ManipulatedVariables(1).Max = num_max_horas_riego * caudal_mm; 
nlobj.ManipulatedVariables(1).MinECR = 0; 
nlobj.ManipulatedVariables(1).MaxECR = 0;

nlobj.OutputVariables(1).Min = R_min; 
nlobj.OutputVariables(1).Max = Rfc;    
nlobj.OutputVariables(1).MinECR = 0.001; 
nlobj.OutputVariables(1).MaxECR = 1; 

nlobj.Weights.OutputVariables = Qm; 
nlobj.Weights.ManipulatedVariables = Rm;

fprintf('Iniciando Análisis de Ruido de Proceso...\n');

for idx = 1:length(niveles_proceso)
    factor_proceso = niveles_proceso(idx);
    fprintf('\n---> Evaluando Ruido de Proceso: ±%d mm/día...\n', factor_proceso);
    
    ruido_sensor = 0.05; 
    error_creciente = linspace(0.05, 0.35, Np); 
    
    Agua_Total_vec = zeros(1, N_simulaciones); 
    Estres_vec = zeros(1, N_simulaciones);
    
    for n = 1:N_simulaciones
        R = zeros(1, N_dias+1); I = zeros(1, N_dias); ET = zeros(1, N_dias); 
        Ke = zeros(1, N_dias); Ks = zeros(1, N_dias); fw = zeros(1, N_dias); 
        Kcmax = zeros(1, N_dias); few = zeros(1, N_dias); fc = zeros(1, N_dias); 
        DP = zeros(1, N_dias); 
        R(1) = 166; last_u = 0;
        
        Pred_Lluvia = zeros(N_dias, Np); Pred_ET0 = zeros(N_dias, Np); 
        Pred_u2 = zeros(N_dias, Np); Pred_hrmin = zeros(N_dias, Np); 
        Pred_Kcb = zeros(N_dias, Np); 
        ruido_constante_kcb = 0.05 * randn(1, Np);  

        for k = 1:N_dias
            lluvia_real = Lluvia_extendida(k:k+Np-1)'; 
            et0_real = ET0_extendida(k:k+Np-1)'; 
            u2_real = u2_extendida(k:k+Np-1)'; 
            hrmin_real = hrmin_extendida(k:k+Np-1)';
            ruido = error_creciente .* randn(1, Np); 
            
            % BLINDAJE 1: Evitar climas negativos
            Pred_Lluvia(k, :) = max(0, lluvia_real .* (1 + ruido)); 
            Pred_ET0(k, :) = max(0, et0_real .* (1 + ruido)); 
            Pred_u2(k, :) = max(0, u2_real .* (1 + ruido)); 
            Pred_hrmin(k, :) = max(0, hrmin_real .* (1 + ruido)); 
            
            for j = 1:Np
                dia_anyo = mod((k+j-1) - 1, 365) + 1;
                if dia_anyo < 60 || dia_anyo >= 330
                    Kcb_tab = 0.4; 
                elseif dia_anyo >= 60 && dia_anyo < 90
                    Kcb_tab = 0.55; 
                elseif dia_anyo >= 90 && dia_anyo < 180
                    Kcb_tab = 0.55 + ((0.65 - 0.55) / 90) * (dia_anyo - 90); 
                else
                    Kcb_tab = 0.65; 
                end
                
                Kcb_base_ruido = Kcb_tab * (1 + ruido_constante_kcb(j));
                ajuste_pred = 0; 
                if dia_anyo >= 180 && dia_anyo < 330
                    ajuste_pred = (0.04*(Pred_u2(k, j)-2) - 0.004*(Pred_hrmin(k, j)-45)) * (h/3)^0.3; 
                end
                Pred_Kcb(k, j) = max(0, Kcb_base_ruido + ajuste_pred);
            end
        end

        for k = 1:N_dias
            % BLINDAJE 2: Limpiar Reserva de basurilla matemática
            R(k) = max(0, real(R(k))); 
            R_medida = max(0, R(k)*(1 + (ruido_sensor * randn()))); 
            
            perturbaciones_futuras = [Pred_Lluvia(k,:)', Pred_ET0(k,:)', Pred_u2(k,:)', Pred_hrmin(k,:)', Pred_Kcb(k,:)'];
            
            [u_opt, ~, info] = nlmpcmove(nlobj, R_medida, last_u, R_objetivo*ones(Np,1), perturbaciones_futuras); 
            
            I(k) = max(0, real(u_opt)); 
            last_u = I(k); 
            
            % BLINDAJE 3: Proteger el cálculo de "fc" de números imaginarios
            Kcb_real_seguro = max(Kcb(k), kcmin + 0.01); 
            Kcmax(k) = max(1.2 + (0.04*(u2(k)-2) - 0.004*(hrmin(k)-45)) * (h/3)^0.3, Kcb_real_seguro+0.05); 
            
            base_fc = (Kcb_real_seguro - kcmin) / (Kcmax(k) - kcmin);
            fc(k) = max(0, min(0.99, base_fc ^ (1 + 0.5*h))); 
            
            if Precipitacion(k) > 3.5
                fw(k) = 1.0; 
            elseif I(k) > 0
                fw(k) = 0.3*(1-((2/3)*fc(k))); 
            else
                if k == 1, fw(k) = 0.3*(1-((2/3)*fc(k))); 
                else, fw(k) = fw(k-1); end
            end
            
            if R(k) >= R_min
                Ks(k) = 1; 
            else
                Ks(k) = max(0, R(k) / R_min); 
            end
            
            Ke(k) = min(1-fc(k), fw(k)) * Kcmax(k); 
            ET(k) = (Ks(k)*Kcb(k) + Ke(k)) * et0(k); 
            DP(k) = do * ((exp(a * R(k)) - 1)/(exp(a * Rfc) - 1));
            
            % Aplicar el ruido brutal de proceso con seguridad
            R(k+1) = max(0, min(Rmax, R(k) + I(k) + Precipitacion(k) - ET(k) - DP(k) + (factor_proceso * randn())));
        end
        Agua_Total_vec(n) = sum(I); 
        Estres_vec(n) = sum(R(1:N_dias) < R_min);
    end
    Agua_Media(idx) = mean(Agua_Total_vec); 
    Estres_Medio(idx) = mean(Estres_vec);
end

% =========================================================================
% GRÁFICAS FINALES (FORMATO GIGANTE)
% =========================================================================
figure('Name', 'Sensibilidad Proceso', 'Color', 'w', 'Position', [100, 100, 1100, 500]);

% --- Subplot 1: Impacto en Consumo de Agua ---
subplot(1,2,1); 
plot(niveles_proceso, Agua_Media, '-o', 'LineWidth', 3, 'MarkerSize', 10, 'MarkerFaceColor', '#0072BD', 'Color', '#0072BD');
xlabel('Ruido de proceso diario (mm)', 'FontSize', 13, 'FontWeight', 'bold'); 
ylabel('Agua Total Riego (mm)', 'FontSize', 13, 'FontWeight', 'bold'); 
title('Impacto en Consumo', 'FontSize', 15, 'FontWeight', 'bold'); 
grid on; 
xticks(niveles_proceso);
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

% --- Subplot 2: Impacto en Seguridad del Cultivo ---
subplot(1,2,2); 
plot(niveles_proceso, Estres_Medio, '-s', 'LineWidth', 3, 'MarkerSize', 10, 'MarkerFaceColor', '#D95319', 'Color', '#D95319');
xlabel('Ruido de proceso diario (mm)', 'FontSize', 13, 'FontWeight', 'bold'); 
ylabel('Días bajo estrés', 'FontSize', 13, 'FontWeight', 'bold'); 
title('Impacto en Seguridad', 'FontSize', 15, 'FontWeight', 'bold'); 
grid on; 
xticks(niveles_proceso); 
ylim([-0.5, max(max(Estres_Medio)+1, 3)]);
set(gca, 'FontSize', 12, 'FontWeight', 'bold')