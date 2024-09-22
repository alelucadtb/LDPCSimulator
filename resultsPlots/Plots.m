dati = readmatrix("fastLDPC2PAM.txt");

X = dati(:, 1);  % Prima colonna (valori X)
Y = dati(:, 2);  % Seconda colonna (valori Y)

% Creare il grafico lineare
plot(Y, X, '-o');  % Linee con punti marcati
xlabel('E_b/N_0', 'FontName', 'Times New Roman');
ylabel('BER', 'FontName', 'Times New Roman');
title('Typical performance of an LDPC code with 8-PAM', 'FontName', 'Times New Roman','FontSize', 16);
grid on;  % Aggiungere la griglia