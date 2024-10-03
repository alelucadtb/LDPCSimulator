dati1 = readmatrix("fastComplexity.txt");
dati2 = readmatrix("slowComplexity.txt");

fastOperations = dati1(:, 2);  % Prima colonna (valori X)
fastSNR = dati1(:, 1);  % Seconda colonna (valori Y)
slowOperations = dati2(:, 2);  % Prima colonna (valori X)
slowSNR = dati2(:, 1);  % Seconda colonna (valori Y)

fastSNR_dB = 10 * log10(fastSNR);
slowSNR_dB = 10 * log10(slowSNR);


% Creare il grafico lineare
plot(fastSNR_dB, fastOperations, '-r', 'LineWidth', 2);  % Linee con punti marcati
hold on
plot(slowSNR_dB, slowOperations, '--b', 'LineWidth', 2);
xlabel('E_b/N_0 [dB]', 'FontName', 'Times New Roman');
ylabel('Sum of costs for the operations', 'FontName', 'Times New Roman');
title('Slow vs Fast decoding cycle LDPC code with 8-PAM', 'FontName', 'Times New Roman','FontSize', 16);
legend('Fast Cycle', 'Slow Cycle');
grid on;  % Aggiungere la griglia