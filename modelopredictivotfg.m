% Modelo: Reserva vs Riego, precipitación, evapotranspiración y drenaje (Cultivo de olivos adultos + Riego predictivo por Goteo en Osuna (Sevilla) , lluvia, et0... del 2025 sacada de la estacion de Osuna)
clear, clc; close all;

% NÚMERO DE SIMULACIONES A EJECUTAR
N_simulaciones = 1;

% 0. Leer el archivo Excel
datos = readtable('Osuna3.csv');
% 0.1 Extraer la columna de precipitación
Precipitacion = flipud(datos.Se11Precip);
% 0.2 Extraer la columna de et0
et0 = flipud(datos.Se11ETo);
% 0.3 Extraer la columna de u2
u2 = flipud(datos.Se11VelViento);
% 0.4 Extraer la columna de hrmin
hrmin = flipud(datos.Se11HumMin);

% 1. Parámetros de la simulación (1 año)
N_dias = 365;             


% 2. Parámetros del tipo de terreno (Franco Arcillo Limoso)
theta_fc = 0.34;           % Contenido de humedad en el suelo a capacidad de campo cuadro 19 pg 167 pdf
theta_wp = 0.20;           % Contenido de humedad en el punto de marchitez permanente cuadro 19 pg 167 pdf
theta_sat = 0.45;          % Contenido de agua a saturación (lo vi en internet)

% 3. Parámetros de la Planta (Olivos Adultos)
Zr = 1.20;                 % Profundidad de las raíces cuadro 22 pg 186 pdf ver nota 1
h=3; %Altura media del árbol en m
kcmin=0.15; %valor mínimo de Kc para suelo sin cobertura y seco [≈ 0,15 - 0,20] pg 172 manual FAO

% 4. Parámetros del cultivo
num_goteros = 4; % por arbol
caudal_got = 4;% %L/h
caudal_total = num_goteros * caudal_got; % L/h por arbol
sup_olivo = 12; %metros cuadrados marco_plantacion (8*5) * fw_goteo (0.3) daria 12 
caudal_mm = caudal_total/sup_olivo; % L/h en un m2 o %% mm/h mejor%% %AEMET 1mm = 1 L/m^2
num_max_horas_riego = 8; 

% 5. Cálculo de la reserva máxima, de la reserva a capacidad de campo y de la R_min 
Rfc = 1000 * (theta_fc - theta_wp) * Zr; % Da 168 mm
Rmax = 1000 * (theta_sat - theta_wp) * Zr; % Da 300 mm
p=0.65; %fraccion de agotamiento
R_min = Rfc * (1 - p); %limite inferior de la banda donde se quiere situar la humedad
R_objetivo = (Rfc + R_min) / 2; % Objetivo mitad de banda

% 6. Parámetros de percolación
do = 1; % Drenaje en mm/día cuando el suelo está a Capacidad de Campo (Rfc)
a = 0.03; % Factor de forma de la curva exponencial



% 8. Parámetros del controlador
Qm = 0.1; % peso para el error de referencia. Hacerla pequeña para priorizar mantenernos en la banda
Rm = 1; % peso para el coste de control

% 9. Ruido del sensor de humedad
ruido_sensor = 0.05; % ruido en la medida

% ===============================================================
% 10. MATRIZ DE PREDICCIÓN Y KCB BASE 
% ===============================================================
Np = 7; % Horizonte de predicción (7 días)
Nc = 7; % Horizonte de control <= Np

% 10.2 Añadimos 7 días extra al final copiando la última semana (para que el bucle no falle en diciembre)
Lluvia_extendida = [Precipitacion; Precipitacion(end-Np+1 : end)];
ET0_extendida    = [et0; et0(end-Np+1 : end)];
u2_extendida     = [u2; u2(end-Np+1 : end)];
hrmin_extendida  = [hrmin; hrmin(end-Np+1 : end)];

%errores
error_creciente = linspace(0.05, 0.35, Np); % Vector de error creciente: del 5% (0.05) el día 1, al 35% (0.35) el día 7 para las predicciones climáticas (lluvia, et0, u2, hrmin)


% 10.3 CÁLCULO PREVIO DEL Kcb REAL (Para todo el año)
Kcb = zeros(1, N_dias);    % Vector de kcb
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

% 10.4 Extender Kcb para que el horizonte predictivo no falle al final del año y asumimos que los últimos días del año y los primeros del siguiente mantienen la tabla base (0.4)
Kcb_extendida = [Kcb, 0.4 * ones(1, Np)];

% ===============================================================
% 11. CONFIGURACIÓN PROFESIONAL DEL CONTROLADOR NLMPC
% ===============================================================
nx = 1; % 1 Estado: Reserva (R)
ny = 1; % 1 Salida: Reserva (R)

% 11.1 Le decimos a MATLAB que en el vector de entradas,el Riego (MV) va en la posición 1, y el clima (MD) va en las posiciones del 2 al 6.
nlobj = nlmpc(nx, ny, 'MV', 1, 'MD', [2 3 4 5 6]);

% 11.2 Tiempos y Horizontes
nlobj.Ts = 1;                % Control diario
nlobj.PredictionHorizon = Np; % Mira Np días adelante (Np)
nlobj.ControlHorizon = Nc;    % Planifica Nc movimientos de riego (Nc)

% 11.3 Vinculamos tu archivo olivo_model.m
nlobj.Model.StateFcn = "olivo_model";
nlobj.Model.IsContinuousTime = false; % Modelo en tiempo discreto

% 11.4 RESTRICCIONES
% 11.4.1 Riego: Restricción DURA (La capacidad del riego)
nlobj.ManipulatedVariables(1).Min = 0; 
nlobj.ManipulatedVariables(1).Max = num_max_horas_riego * caudal_mm; 
nlobj.ManipulatedVariables(1).MinECR = 0; % 0 = Inquebrantable
nlobj.ManipulatedVariables(1).MaxECR = 0;

% 11.4.2 Reserva: Restricción BLANDA. En MATLAB, el peso de una restricción blanda es el INVERSO de su parámetro ECR.
nlobj.OutputVariables(1).Min = R_min;  % No bajar de 58.8 mm
nlobj.OutputVariables(1).Max = Rfc;    % No pasar de 168 mm
nlobj.OutputVariables(1).MinECR = 1;   % multa grande por pasarnos por abajo de la R_min
nlobj.OutputVariables(1).MaxECR = 1; %si nos pasamos por culpa de las precipitaciones no es grave

% 11.5 PESOS (Prioridades del controlador)
% Queremos que el nivel sea prioritario (Qm), pero que ahorre agua (Rm)
nlobj.Weights.OutputVariables = Qm; 
nlobj.Weights.ManipulatedVariables = Rm;

% =========================================================================
% INICIO DEL BUCLE DE SIMULACIONES MÚLTIPLES
% =========================================================================
Agua_Total_Riego_vec = zeros(1, N_simulaciones);
Agua_Total_Drenada_vec = zeros(1, N_simulaciones);
Dias_En_Estres_vec = zeros(1, N_simulaciones);
Dias_Arriba_Rfc_vec = zeros(1, N_simulaciones);

fprintf('Ejecutando %d simulaciones del CONTROLADOR PREDICTIVO MPC...\n', N_simulaciones);

for n = 1:N_simulaciones
    % Reinicio de variables para cada simulacion
    R = zeros(1, N_dias+1);    % Vector de la reserva de agua
    I = zeros(1, N_dias);      % Vector de riego efectivo
    ET = zeros(1, N_dias);     % Vector de ET
    Ke = zeros(1, N_dias);    % Vector de ke
    Ks = zeros(1, N_dias);    % Vector de ks
    fw = zeros(1, N_dias);    % Vector de fw
    Kcmax = zeros(1, N_dias);    % Vector de kc,max
    few = zeros(1, N_dias);    % Vector de few
    fc = zeros(1, N_dias);    % Vector de fc
    DP = zeros(1, N_dias);    %Vector de DP
    R(1) = 166;                 % Suponemos que empezamos el año con 166mm
    % Variables auxiliares para el bucle
    last_u = 0;

    % Generar Matrices de Predicción con ruidos distintos para ESTA simulación
    Pred_Lluvia = zeros(N_dias, Np);
    Pred_ET0    = zeros(N_dias, Np);
    Pred_u2     = zeros(N_dias, Np);
    Pred_hrmin  = zeros(N_dias, Np);
    Pred_Kcb    = zeros(N_dias, Np);
    ruido_constante_kcb = 0.05 * randn(1, Np);  %Error CONSTANTE (ej. 5%) para la tabla del Kcb base

% 10.5 Extraemos los datos "perfectos"
for k = 1:N_dias
    lluvia_real = Lluvia_extendida(k : k + Np - 1)';
    et0_real    = ET0_extendida(k : k + Np - 1)';
    u2_real     = u2_extendida(k : k + Np - 1)';
    hrmin_real  = hrmin_extendida(k : k + Np - 1)';
    
    % 10.5.1 Generamos un ruido aleatorio entre -1 y 1 (randn) multiplicado por el error creciente
    ruido = error_creciente .* randn(1, Np); 
    
    % 10.5.2 Aplicamos el ruido a las predicciones (sumando el porcentaje de error)
    Pred_Lluvia(k, :) = lluvia_real .* (1 + ruido);
    Pred_ET0(k, :)    = et0_real .* (1 + ruido);
    Pred_u2(k, :)     = u2_real .* (1 + ruido);
    Pred_hrmin(k, :)  = hrmin_real .* (1 + ruido);


    % 10.5.3 CÁLCULO DEL Pred_Kcb CON RUIDOS
    for j = 1:Np
        dia_futuro = k + j - 1;
        dia_anyo = mod(dia_futuro - 1, 365) + 1;
        
        % Kcb base de tabla
        if dia_anyo < 60 || dia_anyo >= 330
            Kcb_tab = 0.4;
        elseif dia_anyo >= 60 && dia_anyo < 90
            Kcb_tab = 0.55; 
        elseif dia_anyo >= 90 && dia_anyo < 180
            Kcb_tab = 0.55 + ((0.65 - 0.55) / 90) * (dia_anyo - 90); 
        else 
            Kcb_tab = 0.65; 
        end
        
        % Ruido constante sobre el valor de tabla
        Kcb_base_ruido = Kcb_tab * (1 + ruido_constante_kcb(j));
        
        % Ajuste climático usando las predicciones con ruido creciente
        ajuste_pred = 0;
        if dia_anyo >= 180 && dia_anyo < 330
            ajuste_pred = (0.04*(Pred_u2(k, j)-2) - 0.004*(Pred_hrmin(k, j)-45)) * (h/3)^0.3;
        end
        
        Pred_Kcb(k, j) = max(0, Kcb_base_ruido + ajuste_pred);
    end
end


% =========================================================================
% 12. MODELO: SIMULACIÓN DEL OLIVAR CON CONTROL PREDICTIVO (MPC)
% =========================================================================
for k = 1:N_dias
    % ---------------------------------------------------------------------
    % A. EL CEREBRO DEL MPC DECIDE EL RIEGO DE HOY
    % ---------------------------------------------------------------------
    % Simulación del sensor con ruido
    % El sensor lee la reserva real pero con un pequeño error aleatorio
    R_medida = R(k)*(1 + (ruido_sensor * randn()));
    R_medida = max(0, R_medida); % Por si el ruido hace que el sensor lea negativo
    % 1. Empaquetamos las predicciones de los próximos 7 días (7 filas x 5 columnas)
    % Columnas: [Lluvia, ET0, u2, hrmin, kcb]
    perturbaciones_futuras = [Pred_Lluvia(k,:)', Pred_ET0(k,:)', Pred_u2(k,:)', Pred_hrmin(k,:)', Pred_Kcb(k,:)'];
    
    % 2. El optimizador calcula el mejor movimiento
    % Le pasamos: el objeto nlobj, el nivel actual R(k), el riego de ayer (last_u), el objetivo a alcanzar y las previsiones del clima.
    [u_opt, ~, info] = nlmpcmove(nlobj, R_medida, last_u, R_objetivo*ones(Np,1), perturbaciones_futuras); 
    
    % 3. Aplicamos el riego decidido (asegurando que no sea negativo por fallos numéricos)
    I(k) = max(0, u_opt);
    last_u = I(k); % Guardamos lo que hemos regado para que el MPC lo sepa mañana
  
    % 4. Aplicamos Kcmax y fc
    Kcmax(k) = max(1.2 + (0.04*(u2(k)-2) - 0.004*(hrmin(k)-45)) * (h/3)^0.3, Kcb(k)+0.05); 
    fc(k)= ((Kcb(k) - kcmin) / (Kcmax(k) - kcmin)) ^ (1 + 0.5*h);
    if fc(k) < 0, fc(k) = 0; elseif fc(k) > 0.99, fc(k) = 0.99; end
    
    % 5. Fracción de suelo humedecida (fw)
    umbral_lluvia = 3.5; 
    fwgoteo = 0.3; 
    
    if Precipitacion(k) > umbral_lluvia 
        fw(k) = 1.0;  % Lluvia fuerte
    elseif I(k) > 0
        fw(k) = fwgoteo*(1-((2/3)*fc(k))); % El MPC ha decidido regar hoy
    else 
        if k == 1
            fw(k) = fwgoteo*(1-((2/3)*fc(k))); 
        else
            fw(k) = fw(k-1);  
        end
    end
    
    %Cálculo del coeficiente de estrés (Ks) para el SISTEMA REAL
    if R(k) >= R_min
        Ks(k) = 1; % No hay estrés
    else
        Ks(k) = R(k) / R_min; % Hay estrés, la planta bebe menos
    end

    % 6. Evaporación del suelo y ETc
    few(k) = min(1-fc(k), fw(k)); 
    Ke(k) = few(k) * Kcmax(k);
    ET(k) = (Ks(k)*Kcb(k) + Ke(k)) * et0(k);
    
    % 7. Percolación profunda exponencial
    DP(k) = do * ((exp(a * R(k)) - 1)/(exp(a * Rfc) - 1));
    
    % 8. BALANCE DE MASAS DIARIO
    ruido_proceso = 1.0 * randn();
    R(k+1) = R(k) + I(k) + Precipitacion(k) - ET(k) - DP(k) + ruido_proceso;
    
    % 9. Limites de saturación física del suelo
    if R(k+1) < 0
        R(k+1) = 0; 
    elseif R(k+1) > Rmax 
        R(k+1) = Rmax; 
    end
end

    Agua_Total_Riego_vec(n) = sum(I);
    Agua_Total_Drenada_vec(n) = sum(DP);
    Dias_En_Estres_vec(n) = sum(R(1:N_dias) < R_min);
    Dias_Arriba_Rfc_vec(n) = sum(R(1:N_dias) > Rfc);

end


% =========================================================================
% CÁLCULO DE MEDIAS Y DESVIACIONES ESTÁNDAR
% =========================================================================
fprintf('\n=== RESULTADOS MEDIOS TRAS %d SIMULACIONES: MPC ===\n', N_simulaciones);
fprintf('Agua total aplicada en riego: %.2f ± %.2f mm\n', mean(Agua_Total_Riego_vec), std(Agua_Total_Riego_vec));
fprintf('Agua total perdida por drenaje: %.2f ± %.2f mm\n', mean(Agua_Total_Drenada_vec), std(Agua_Total_Drenada_vec));
fprintf('Número de días bajo el umbral de estrés: %.1f ± %.1f días\n', mean(Dias_En_Estres_vec), std(Dias_En_Estres_vec));
fprintf('Número de días arriba de capacidad de campo: %.1f ± %.1f días\n', mean(Dias_Arriba_Rfc_vec), std(Dias_Arriba_Rfc_vec));

% =========================================================================
% PLOTEO DE RESULTADOS GRÁFICOS (Muestra la última simulación)
% =========================================================================

% --- GRÁFICA 1: RESERVA DE AGUA ---
figure('Name', 'Reserva de Agua', 'Position', [100, 100, 900, 500]);
plot(1:N_dias+1, R, 'b', 'LineWidth', 2.5); hold on;
yline(Rmax, 'r--', 'Reserva Máx. (300 mm)', 'LineWidth', 2, 'LabelHorizontalAlignment', 'left', 'FontSize', 11, 'FontWeight', 'bold');
yline(Rfc, 'g--', 'Capacidad de campo (168 mm)', 'LineWidth', 2, 'LabelHorizontalAlignment', 'right', 'FontSize', 11, 'FontWeight', 'bold');
yline(R_objetivo, 'k--', 'Objetivo de Riego', 'LineWidth', 2, 'LabelHorizontalAlignment', 'left', 'FontSize', 11, 'FontWeight', 'bold');
yline(R_min, 'm--', 'Reserva Mínima (58.8 mm)', 'LineWidth', 2, 'LabelHorizontalAlignment', 'right', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Reserva (mm)', 'FontSize', 14, 'FontWeight', 'bold'); 
xlabel('Días (k)', 'FontSize', 14, 'FontWeight', 'bold'); 
ylim([0, Rmax + 50]);
title('Evolución de la Reserva (Control MPC)', 'FontSize', 16, 'FontWeight', 'bold'); 
grid on;
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

% --- GRÁFICA 2: PERCOLACIÓN PROFUNDA ---
figure('Name', 'Percolación Profunda', 'Position', [150, 150, 900, 500]);
plot(R(1:N_dias), DP, 'b.', 'MarkerSize', 12); hold on;
xline(Rfc, 'g--', 'Capacidad de campo', 'LineWidth', 2, 'FontSize', 11, 'FontWeight', 'bold'); 
ylabel('Percolación DP(k) (mm)', 'FontSize', 14, 'FontWeight', 'bold'); 
xlabel('Reserva R(k) (mm)', 'FontSize', 14, 'FontWeight', 'bold');
title('Relación entre Reserva y Drenaje', 'FontSize', 16, 'FontWeight', 'bold'); 
grid on;
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

% --- GRÁFICA 3: ESTRATEGIA DE RIEGO ---
figure('Name', 'Estrategia de Riego MPC', 'Position', [200, 200, 900, 500]);
hold on;
bar(1:N_dias, Precipitacion, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', 'DisplayName', 'Lluvia');
bar(1:N_dias, I, 'FaceColor', [0 0.4 0.8], 'EdgeColor', 'none', 'DisplayName', 'Riego Optimo MPC');
ylabel('Volumen (mm/día)', 'FontSize', 14, 'FontWeight', 'bold'); 
xlabel('Días (k)', 'FontSize', 14, 'FontWeight', 'bold');
title('Estrategia de Riego Inteligente: Controlador Predictivo', 'FontSize', 16, 'FontWeight', 'bold');
lgd = legend('Location', 'northeast'); 
lgd.FontSize = 12; 
lgd.FontWeight = 'bold';
xlim([1, N_dias]); 
grid on;
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
