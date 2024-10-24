dati1 = readmatrix("8PAMInterleavingDepth.txt");
dati2 = readmatrix("8PAMInterleavingDepth_v2.txt");
dati3 = readmatrix("8PAMInterleavingDepth_v3.txt")

Depth1 = dati1(:, 2);  % Prima colonna (valori X)
BER1 = dati1(:, 1);  % Seconda colonna (valori Y)
Depth2 = dati2(:, 2);
BER2 = dati2(:, 1);
Depth3 = dati3(:, 2);
BER3 = dati3(:, 1);

% Creare il grafico lineare
figure;
loglog(Depth1, BER1, '-r', 'LineWidth', 2);
hold on
loglog(Depth2, BER2, '--b', 'LineWidth', 2);
hold on
loglog(Depth3, BER3, '--g', 'LineWidth', 2);
xlabel('Interleaving Depth', 'FontName', 'Times New Roman');
ylabel('BER', 'FontName', 'Times New Roman');
title('Different interleaving approach with 8-PAM', 'FontName', 'Times New Roman','FontSize', 16);
legend('Interleaving\_v1', 'Interleaving\_v2', 'Interleaving\_v3');
grid on;  % Aggiungere la griglia