%Control riego del técnico 3 niveles
%Olivar adulto en Osuna, simulacion con los datos del 2025 y simulando N_simulaciones
clear, clc; close all;

% SELECCIÓN DEL NIVEL DEL TÉCNICO AGRÓNOMO 
% 1 = Nivel Básico (Solo Reserva R_medida)
% 2 = Nivel Intermedio (Reserva R_medida + Predicción de Lluvia)
% 3 = Nivel Avanzado (Reserva R_medida + Predicción de Lluvia + Cálculo ETc)
nivel_tecnico = 3;

% NÚMERO DE SIMULACIONES A EJECUTAR
N_simulaciones = 100;

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

R_objetivo = (Rfc + R_min) / 2;

% 6. Parámetros de percolación
do = 1; % Drenaje en mm/día cuando el suelo está a Capacidad de Campo (Rfc)
a = 0.03; % Factor de forma de la curva exponencial

% 7. Condiciones Iniciales
ruido_sensor = 0.05; % 5% de ruido en la medida del sensor

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

% Vectores para almacenar los resultados de las N simulaciones
Agua_Total_Riego_vec = zeros(1, N_simulaciones);
Agua_Total_Drenada_vec = zeros(1, N_simulaciones);
Dias_En_Estres_vec = zeros(1, N_simulaciones);
Dias_Arriba_Rfc_vec = zeros(1, N_simulaciones);

fprintf('Ejecutando %d simulaciones del Nivel Técnico %d...\n', N_simulaciones, nivel_tecnico);

% =========================================================================
% INICIO DEL BUCLE DE SIMULACIONES MÚLTIPLES
% =========================================================================
for n = 1:N_simulaciones
    
    % IMPORTANTE: Reiniciamos las variables CADA simulacion para que no se acumulen
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
    R(1) = 50;                 % Suponemos que empezamos el año con 166mm

% =========================================================================
%  BUCLE PRINCIPAL DE SIMULACIÓN DEL ENTORNO
% =========================================================================
for k = 1:N_dias
    R_medida = R(k) * (1 + (ruido_sensor * randn()));
    R_medida = max(0, R_medida); % Evitar lecturas físicas negativas

    ruido_clima = 0.05 * randn();

    lluvia_esperada_hoy = max(0, Precipitacion(k) * (1 + ruido_clima));
    et0_estimado        = max(0, et0(k) * (1 + ruido_clima));
    u2_estimado         = max(0, u2(k) * (1 + ruido_clima));
    hrmin_estimado      = max(0, hrmin(k) * (1 + ruido_clima));

    ruido_kcb = 0.05 * randn();

    switch nivel_tecnico
        case 1
            % NIVEL 1: TÉCNICO BÁSICO 
            if R_medida < R_objetivo
                I(k) = R_objetivo - R_medida; 
            else
                I(k) = 0;
            end
            
        case 2
            % NIVEL 2: TÉCNICO INTERMEDIO 
            if R_medida < R_objetivo
                I(k) = R_objetivo - R_medida - lluvia_esperada_hoy;
            else
                I(k) = 0;
            end
            
        case 3
            % NIVEL 3: TÉCNICO AVANZADO
            ajuste_estimado = 0;
            if k < 60 || k >= 330
                Kcb_tab = 0.4;
            elseif k >= 60 && k < 90
                Kcb_tab = 0.55; 
            elseif k >= 90 && k < 180
                Kcb_tab = 0.55 + ((0.65 - 0.55) / 90) * (k - 90); 
            else 
                Kcb_tab = 0.65; 
                ajuste_estimado = (0.04*(u2_estimado-2) - 0.004*(hrmin_estimado-45)) * (h/3)^0.3;
            end
                Kcb_base_ruido = Kcb_tab * (1 + ruido_kcb);
                Kcb_estimado = max(0, Kcb_base_ruido + ajuste_estimado);

            Kcmax_teorica = max(1.2 + (0.04*(u2_estimado-2) - 0.004*(hrmin_estimado-45)) * (h/3)^0.3, Kcb_estimado+0.05);
            fc_teorica = ((Kcb_estimado - kcmin) / (Kcmax_teorica - kcmin)) ^ (1 + 0.5*h);
            if fc_teorica < 0, fc_teorica = 0; elseif fc_teorica > 0.99, fc_teorica = 0.99; end

            if k == 1, fw_estimada = 0.3; else, fw_estimada = fw(k-1); end
            few_teorica = min(1-fc_teorica, fw_estimada);
            Ke_teorica = few_teorica * Kcmax_teorica;
             
            ET_estimada_tecnico = (Kcb_estimado + Ke_teorica) * et0_estimado;

            if R_medida < R_objetivo
                I(k) = R_objetivo - R_medida + ET_estimada_tecnico - lluvia_esperada_hoy;
            else
                I(k) = 0;
            end
    end

    I_maximo=caudal_mm*num_max_horas_riego;

    if I(k) > I_maximo
            I(k) = I_maximo;
    end
    I(k) = max(0, I(k));

    % Aplicamos Kcmax y fc
    Kcmax(k) = max(1.2 + (0.04*(u2(k)-2) - 0.004*(hrmin(k)-45)) * (h/3)^0.3, Kcb(k)+0.05); 
    fc(k)= ((Kcb(k) - kcmin) / (Kcmax(k) - kcmin)) ^ (1 + 0.5*h);
    if fc(k) < 0, fc(k) = 0; elseif fc(k) > 0.99, fc(k) = 0.99; end

    %Fracción de suelo humedecida (fw)
    umbral_lluvia = 3.5; 
    fwgoteo = 0.3; 

    if Precipitacion(k) > umbral_lluvia 
        fw(k) = 1.0;  % Si llueve significativamente (se riegue o no) se moja todo el campo
    elseif I(k) > 0
        fw(k) = fwgoteo*(1-((2/3)*fc(k))); % El MPC ha decidido regar hoy
    else % No hay ni riego ni lluvia significativa
        if k == 1
            fw(k) = fwgoteo*(1-((2/3)*fc(k))); % El MPC ha decidido regar hoy
        else
            fw(k) = fw(k-1);  % Mantiene el valor del día anterior
        end
    end


    %Cálculo del coeficiente de estrés (Ks) para el SISTEMA REAL
    if R(k) >= R_min
        Ks(k) = 1; % No hay estrés
    else
        Ks(k) = R(k) / R_min; % Hay estrés, la planta bebe menos
    end

    %  Evaporación del suelo y ETc
    few(k) = min(1-fc(k), fw(k)); 
    Ke(k) = few(k) * Kcmax(k);
    ET(k) = (Ks(k)*Kcb(k) + Ke(k)) * et0(k);

     %  Percolación profunda exponencial
    DP(k) = do * ((exp(a * R(k)) - 1)/(exp(a * Rfc) - 1));

    ruido_proceso = 1.0 * randn();
    R(k+1) = R(k) + I(k) + Precipitacion(k) - ET(k) - DP(k) + ruido_proceso;
    
    % Límites físicos de saturación del suelo real
    if R(k+1) < 0
        R(k+1) = 0; 
    elseif R(k+1) > Rmax 
        R(k+1) = Rmax; 
    end
end

% Guardar los resultados de esta simulación 'n'
    Agua_Total_Riego_vec(n) = sum(I);
    Agua_Total_Drenada_vec(n) = sum(DP);
    Dias_En_Estres_vec(n) = sum(R(1:N_dias) < R_min);
    Dias_Arriba_Rfc_vec(n) = sum(R(1:N_dias) > Rfc);
    
end 

% =========================================================================
% CÁLCULO DE MEDIAS Y DESVIACIONES ESTÁNDAR
% =========================================================================
nombres_niveles = {'TÉCNICO BÁSICO', 'TÉCNICO INTERMEDIO', 'TÉCNICO AVANZADO'};
nombre_actual = nombres_niveles{nivel_tecnico};

fprintf('\n=== RESULTADOS MEDIOS TRAS %d SIMULACIONES: %s ===\n', N_simulaciones, nombre_actual);
fprintf('Agua total aplicada en riego: %.2f ± %.2f mm\n', mean(Agua_Total_Riego_vec), std(Agua_Total_Riego_vec));
fprintf('Agua total perdida por drenaje: %.2f ± %.2f mm\n', mean(Agua_Total_Drenada_vec), std(Agua_Total_Drenada_vec));
fprintf('Número de días bajo el umbral de estrés: %.1f ± %.1f días\n', mean(Dias_En_Estres_vec), std(Dias_En_Estres_vec));
fprintf('Número de días arriba de capacidad de campo: %.1f ± %.1f días\n', mean(Dias_Arriba_Rfc_vec), std(Dias_Arriba_Rfc_vec));


% =========================================================================
% PLOTEO DE RESULTADOS GRÁFICOS (Muestra la última simulación)
% =========================================================================

% --- GRÁFICA 1: RESERVA DE AGUA ---
figure('Name', ['Reserva de Agua - ' nombre_actual], 'Position', [100, 100, 900, 500]);
plot(1:N_dias+1, R, 'b', 'LineWidth', 2.5); 
hold on;
yline(Rmax, 'r--', 'Reserva Máx. (300 mm)', 'LineWidth', 2, 'LabelHorizontalAlignment', 'left', 'FontSize', 11, 'FontWeight', 'bold');
yline(Rfc, 'g--', 'Capacidad de campo (168 mm)', 'LineWidth', 2, 'LabelHorizontalAlignment', 'right', 'FontSize', 11, 'FontWeight', 'bold');
yline(R_objetivo, 'k--', 'Objetivo de Riego', 'LineWidth', 2, 'LabelHorizontalAlignment', 'left', 'FontSize', 11, 'FontWeight', 'bold');
yline(R_min, 'm--', 'Reserva Mínima (58.8 mm)', 'LineWidth', 2, 'LabelHorizontalAlignment', 'right', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Reserva (mm)', 'FontSize', 14, 'FontWeight', 'bold'); 
xlabel('Días (k)', 'FontSize', 14, 'FontWeight', 'bold'); 
ylim([0, Rmax + 50]);
title(['Evolución de la Reserva (Control Tradicional: ' nombre_actual ')'], 'FontSize', 16, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12, 'FontWeight', 'bold'); % Hace grandes los números de los ejes

% --- GRÁFICA 2: ESTRATEGIA DE RIEGO ---
figure('Name', ['Estrategia de Riego - ' nombre_actual], 'Position', [150, 150, 900, 500]);
hold on;
bar(1:N_dias, Precipitacion, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', 'DisplayName', 'Lluvia');
bar(1:N_dias, I, 'FaceColor', [0.85 0.33 0.1], 'EdgeColor', 'none', 'DisplayName', ['Riego ' nombre_actual]);
ylabel('Volumen (mm/día)', 'FontSize', 14, 'FontWeight', 'bold'); 
xlabel('Días (k)', 'FontSize', 14, 'FontWeight', 'bold');
title(['Estrategia de Riego: ' nombre_actual], 'FontSize', 16, 'FontWeight', 'bold');
lgd = legend('Location', 'northeast');
lgd.FontSize = 12; % Hace grande la letra de la leyenda
lgd.FontWeight = 'bold';
xlim([1, N_dias]); 
grid on;
set(gca, 'FontSize', 12, 'FontWeight', 'bold'); % Hace grandes los números de los ejes