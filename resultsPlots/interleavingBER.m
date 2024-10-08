dati1 = readmatrix("8PAMInterleavingDepth.txt");
dati2 = readmatrix("4PAMInterleavingDepth.txt");

Depth1 = dati1(:, 2);  % Prima colonna (valori X)
BER1 = dati1(:, 1);  % Seconda colonna (valori Y)
Depth2 = dati2(:, 2);
BER2 = dati2(:, 1);

% Creare il grafico lineare
figure;
loglog(Depth1, BER1, '-r', 'LineWidth', 2);
hold on
loglog(Depth2, BER2, '-b', 'LineWidth', 2);
xlabel('Interleaving Depth', 'FontName', 'Times New Roman');
ylabel('BER', 'FontName', 'Times New Roman');
title('Slow vs Fast decoding cycle LDPC code with 8-PAM', 'FontName', 'Times New Roman','FontSize', 16);
legend('8-PAM', '4-PAM');
grid on;  % Aggiungere la griglia