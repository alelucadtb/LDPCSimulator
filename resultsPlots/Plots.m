dati1 = readmatrix("fastCycleLDPC8PAM.txt");
dati2 = readmatrix("slowCycleLDPC8PAM.txt");

fastBER = dati1(:, 1);  % Prima colonna (valori X)
fastSNR = dati1(:, 2);  % Seconda colonna (valori Y)
slowBER = dati2(:, 1);
slowSNR = dati2(:, 2);

fastSNR_dB = 20 * log10(fastSNR);
slowSNR_dB = 20 * log10(slowSNR);

% Creare il grafico lineare
figure;
plot(fastSNR_dB, fastBER, '-r', 'LineWidth', 2);
hold on
plot(slowSNR_dB, slowBER, '--b', 'LineWidth', 2);
xlabel('E_b/N_0 [dB]', 'FontName', 'Times New Roman');
ylabel('BER', 'FontName', 'Times New Roman');
title('Fast vs Slow decoding cycle LDPC code with 8-PAM', 'FontName', 'Times New Roman','FontSize', 16);
legend('Fast cycle', 'Slow cycle');
grid on;  % Aggiungere la griglia
