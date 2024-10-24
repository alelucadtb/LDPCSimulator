dati1 = readmatrix("fCLDPC8PAMR56.txt");
dati2 = readmatrix("sCLDPC8PAMR56.txt");

fastBER = dati1(:, 1);  % Prima colonna (valori X)
fastSNR = dati1(:, 2);  % Seconda colonna (valori Y)
slowBER = dati2(:, 1);
slowSNR = dati2(:, 2);

fastSNR_dB = 10 * log10(fastSNR);
slowSNR_dB = 10 * log10(slowSNR);
fastBER_dB = 10 * log10(fastBER);
slowBER_dB = 10 * log10(slowBER);

% Creare il grafico lineare
figure;
loglog(fastSNR, fastBER, '-r', 'LineWidth', 2);
hold on
loglog(slowSNR, slowBER, '--b', 'LineWidth', 2);
xlabel('E_b/N_0 [dB]', 'FontName', 'Times New Roman');
ylabel('BER', 'FontName', 'Times New Roman');
title('Fast vs Slow decoding cycle LDPC code with 8-PAM', 'FontName', 'Times New Roman','FontSize', 16);
legend('Fast cycle', 'Slow cycle');
grid on;  % Aggiungere la griglia
