dati = readmatrix("fastLDPC8PAM.txt");

BER = dati(:, 1);  % Prima colonna (valori X)
SNR = dati(:, 2);  % Seconda colonna (valori Y)

% Creare il grafico lineare
plot(SNR, BER, '-o');  % Linee con punti marcati
xlabel('E_b/N_0', 'FontName', 'Times New Roman');
ylabel('BER', 'FontName', 'Times New Roman');
title('Fast decoding cycle LDPC code with 8-PAM', 'FontName', 'Times New Roman','FontSize', 16);
grid on;  % Aggiungere la griglia