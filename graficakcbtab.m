% Generar gráfica del Kcb base anual para el TFG
N_dias = 365;
Kcb_base = zeros(1, N_dias);

for k = 1:N_dias
    if k < 60 || k >= 330
        Kcb_base(k) = 0.4;
    elseif k >= 60 && k < 90
        Kcb_base(k) = 0.55; 
    elseif k >= 90 && k < 180
        Kcb_base(k) = 0.55 + ((0.65 - 0.55) / 90) * (k - 90); 
    else 
        Kcb_base(k) = 0.65; 
    end
end

figure('Name', 'Evolución del Kcb base');
plot(1:N_dias, Kcb_base, 'b', 'LineWidth', 2.5);
xlabel('Días del año (k)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Coeficiente basal tabulado (K_{cb,tab})', 'FontSize', 12, 'FontWeight', 'bold');
title('Evolución anual del K_{cb} base para olivar', 'FontSize', 14, 'FontWeight', 'bold');
ylim([0.3 0.75]);
xlim([1 365]);
grid on;
set(gca, 'FontSize', 12);