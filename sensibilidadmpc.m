% =========================================================================
% SCRIPT DE ANÁLISIS DE SENSIBILIDAD AISLADO (SENSOR) - APARTADO 6.4.3
% =========================================================================
clear, clc; close all;

% 1. CONFIGURACIÓN DEL ANÁLISIS DE SENSIBILIDAD
niveles_ruido = [5, 15, 25, 35, 50]; % Porcentajes de ruido a evaluar SOLO PARA EL SENSOR
N_simulaciones = 10; % Simulaciones por cada nivel de ruido (Ajustar según PC)

% Vectores para guardar los resultados finales de las gráficas
Agua_Media_Sensibilidad = zeros(1, length(niveles_ruido));
Estres_Medio_Sensibilidad = zeros(1, length(niveles_ruido));

% 2. LECTURA DE DATOS BASE 
datos = readtable('Osuna3.csv');
Precipitacion = flipud(datos.Se11Precip);
et0 = flipud(datos.Se11ETo);
u2 = flipud(datos.Se11VelViento);
hrmin = flipud(datos.Se11HumMin);

N_dias = 365;             
theta_fc = 0.34;           
theta_wp = 0.20;           
theta_sat = 0.45;          
Zr = 1.20;                 
h=3; 
kcmin=0.15; 
num_goteros = 4; 
caudal_got = 4;
caudal_total = num_goteros * caudal_got; 
sup_olivo = 12; 
caudal_mm = caudal_total/sup_olivo; 
num_max_horas_riego = 8; 

Rfc = 1000 * (theta_fc - theta_wp) * Zr; 
Rmax = 1000 * (theta_sat - theta_wp) * Zr; 
p=0.65; 
R_min = Rfc * (1 - p); 
R_objetivo = (Rfc + R_min) / 2; 

do = 1; 
a = 0.03; 

Qm = 0.1; 
Rm = 1; 

Np = 7; 
Nc = 7; 
Lluvia_extendida = [Precipitacion; Precipitacion(end-Np+1 : end)];
ET0_extendida    = [et0; et0(end-Np+1 : end)];
u2_extendida     = [u2; u2(end-Np+1 : end)];
hrmin_extendida  = [hrmin; hrmin(end-Np+1 : end)];

Kcb = zeros(1, N_dias);    
for k = 1:N_dias
    ajuste = 0;
    if k < 60 || k >= 330
        Kcb_tab = 0.4;
    elseif k >= 60 && k < 90
        Kcb_tab = 0.55; 
    elseif k >= 90 && k < 180
        Kcb_tab = 0.55 + ((0.65 - 0.55) / 90) * (k - 90); 
    else 
        Kcb_tab = 0.65; 
        ajuste = (0.04*(u2(k)-2) - 0.004*(hrmin(k)-45)) * (h/3)^0.3;
    end
    Kcb(k) = Kcb_tab + ajuste;
end
Kcb_extendida = [Kcb, 0.4 * ones(1, Np)];

nx = 1; 
ny = 1; 
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

% =========================================================================
% SÚPER-BUCLE DE SENSIBILIDAD (RECORRE LOS NIVELES DE RUIDO DEL SENSOR)
% =========================================================================
fprintf('Iniciando Análisis de Sensibilidad (Aislado para el Sensor)...\n');

for idx_ruido = 1:length(niveles_ruido)
    ruido_actual = niveles_ruido(idx_ruido) / 100; % Convertir a tanto por uno
    fprintf('\n---> Evaluando escenario con Sensor al %d%% de error...\n', niveles_ruido(idx_ruido));
    
    % --- AISLAMIENTO DE VARIABLES ---
    % 1. ESTE ES EL ÚNICO QUE VARÍA
    ruido_sensor = ruido_actual; 
    
    % 2. ESTOS SE QUEDAN FIJOS EN SU VALOR NOMINAL (5%)
    error_creciente = linspace(0.05, 0.35, Np); 
    factor_proceso = 1.0; 

    Agua_Total_Riego_vec = zeros(1, N_simulaciones);
    Dias_En_Estres_vec = zeros(1, N_simulaciones);
    
    % BUCLE DE MONTECARLO PARA EL RUIDO ACTUAL DEL SENSOR
    for n = 1:N_simulaciones
        R = zeros(1, N_dias+1);    
        I = zeros(1, N_dias);      
        ET = zeros(1, N_dias);     
        Ke = zeros(1, N_dias);    
        Ks = zeros(1, N_dias);    
        fw = zeros(1, N_dias);    
        Kcmax = zeros(1, N_dias);    
        few = zeros(1, N_dias);    
        fc = zeros(1, N_dias);    
        DP = zeros(1, N_dias);    
        R(1) = 166;                 
        last_u = 0;
        
        Pred_Lluvia = zeros(N_dias, Np);
        Pred_ET0    = zeros(N_dias, Np);
        Pred_u2     = zeros(N_dias, Np);
        Pred_hrmin  = zeros(N_dias, Np);
        Pred_Kcb    = zeros(N_dias, Np);
        
        ruido_constante_kcb = 0.05 * randn(1, Np);  % Fijo al 5%

        for k = 1:N_dias
            lluvia_real = Lluvia_extendida(k : k + Np - 1)';
            et0_real    = ET0_extendida(k : k + Np - 1)';
            u2_real     = u2_extendida(k : k + Np - 1)';
            hrmin_real  = hrmin_extendida(k : k + Np - 1)';
            
            ruido = error_creciente .* randn(1, Np); 
            
            Pred_Lluvia(k, :) = lluvia_real .* (1 + ruido);
            Pred_ET0(k, :)    = et0_real .* (1 + ruido);
            Pred_u2(k, :)     = u2_real .* (1 + ruido);
            Pred_hrmin(k, :)  = hrmin_real .* (1 + ruido);
            
            for j = 1:Np
                dia_futuro = k + j - 1;
                dia_anyo = mod(dia_futuro - 1, 365) + 1;
                
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
            % AQUI ENTRA EN JUEGO EL SENSOR VARIABLE
            R_medida = R(k)*(1 + (ruido_sensor * randn()));
            R_medida = max(0, R_medida); 
            
            perturbaciones_futuras = [Pred_Lluvia(k,:)', Pred_ET0(k,:)', Pred_u2(k,:)', Pred_hrmin(k,:)', Pred_Kcb(k,:)'];
            
            [u_opt, ~, info] = nlmpcmove(nlobj, R_medida, last_u, R_objetivo*ones(Np,1), perturbaciones_futuras); 
            
            I(k) = max(0, u_opt);
            last_u = I(k); 
          
            Kcmax(k) = max(1.2 + (0.04*(u2(k)-2) - 0.004*(hrmin(k)-45)) * (h/3)^0.3, Kcb(k)+0.05); 
            fc(k)= ((Kcb(k) - kcmin) / (Kcmax(k) - kcmin)) ^ (1 + 0.5*h);
            if fc(k) < 0, fc(k) = 0; elseif fc(k) > 0.99, fc(k) = 0.99; end
            
            umbral_lluvia = 3.5; 
            fwgoteo = 0.3; 
            
            if Precipitacion(k) > umbral_lluvia 
                fw(k) = 1.0;  
            elseif I(k) > 0
                fw(k) = fwgoteo*(1-((2/3)*fc(k))); 
            else 
                if k == 1
                    fw(k) = fwgoteo*(1-((2/3)*fc(k))); 
                else
                    fw(k) = fw(k-1);  
                end
            end
            
            if R(k) >= R_min
                Ks(k) = 1; 
            else
                Ks(k) = R(k) / R_min; 
            end
            
            few(k) = min(1-fc(k), fw(k)); 
            Ke(k) = few(k) * Kcmax(k);
            ET(k) = (Ks(k)*Kcb(k) + Ke(k)) * et0(k);
            
            DP(k) = do * ((exp(a * R(k)) - 1)/(exp(a * Rfc) - 1));
            
            ruido_proceso = factor_proceso * randn();
            R(k+1) = R(k) + I(k) + Precipitacion(k) - ET(k) - DP(k) + ruido_proceso;
            
            if R(k+1) < 0
                R(k+1) = 0; 
            elseif R(k+1) > Rmax 
                R(k+1) = Rmax; 
            end
        end
        
        Agua_Total_Riego_vec(n) = sum(I);
        Dias_En_Estres_vec(n) = sum(R(1:N_dias) < R_min);
    end
    
    Agua_Media_Sensibilidad(idx_ruido) = mean(Agua_Total_Riego_vec);
    Estres_Medio_Sensibilidad(idx_ruido) = mean(Dias_En_Estres_vec);
end

fprintf('\n=== ANÁLISIS COMPLETADO ===\nGenerando gráficas...\n');

% =========================================================================
% GENERACIÓN DE GRÁFICAS PARA EL APARTADO 6.4.3
% =========================================================================

% 1. Gráfica de Degradación del Consumo de Agua
figure('Name', 'Sensibilidad: Agua Total', 'Color', 'w', 'Position', [100, 100, 850, 500]);
plot(niveles_ruido, Agua_Media_Sensibilidad, '-o', 'LineWidth', 3, 'MarkerSize', 10, 'MarkerFaceColor', '#0072BD', 'Color', '#0072BD');
title('Sensibilidad del NMPC frente al Ruido del Sensor de Humedad', 'FontSize', 15, 'FontWeight', 'bold');
xlabel('Incertidumbre del Sensor (%)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Agua Total Aplicada (mm)', 'FontSize', 13, 'FontWeight', 'bold');
grid on;
xticks(niveles_ruido); 
ylim([min(Agua_Media_Sensibilidad)-10, max(Agua_Media_Sensibilidad)+20]);
set(gca, 'FontSize', 12, 'FontWeight', 'bold'); % Agranda números de los ejes y los pone en negrita

% 2. Gráfica de Robustez (Días de Estrés)
figure('Name', 'Sensibilidad: Seguridad del Cultivo', 'Color', 'w', 'Position', [150, 150, 850, 500]);
plot(niveles_ruido, Estres_Medio_Sensibilidad, '-s', 'LineWidth', 3, 'MarkerSize', 10, 'MarkerFaceColor', '#D95319', 'Color', '#D95319');
title('Robustez del NMPC frente al Ruido del Sensor de Humedad', 'FontSize', 15, 'FontWeight', 'bold');
xlabel('Incertidumbre del Sensor (%)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Días bajo estrés hídrico', 'FontSize', 13, 'FontWeight', 'bold');
grid on;
xticks(niveles_ruido);
ylim([-0.5, max(max(Estres_Medio_Sensibilidad)+1, 3)]);
set(gca, 'FontSize', 12, 'FontWeight', 'bold'); % Agranda números de los ejes y los pone en negrita