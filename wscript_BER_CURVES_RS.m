%   SIMCOM DVB-S Simulator
%   2014/2015 Juan Pablo Cuadro and Loic Veillard

close all
clear all

N_SYMBOLS = 4 * 5 * 1.504 * 1e4;
LABEL_SCENARIOS = {'Theoretical',...
    'Simulated - no FEC',...
    'Simulated - Conv Hard Deco - R = 1/2',...
    'Simulated - Conv Soft Deco - R = 1/2',...
    'Simulated - Conv Punctured - R = 2/3',...
    'Simulated - No Interleave  - RS',...
    'Simulated - Interleaved    - RS'    
    };
N_SCENARIOS     = numel(LABEL_SCENARIOS);

EbN0_dB = -4:1:7;

BER = zeros(N_SCENARIOS, length(EbN0_dB));

% Construct Simulator Class
obj = DVBS_Simulator('nofec', N_SYMBOLS);

% Theoretical BER
BER(1,:) = 0.5*erfc(sqrt(10.^(EbN0_dB/10)));


% profile on 

tic
progressbar(0);

for i = 1:length(EbN0_dB)
    
    % **** Scenarios! **** %
    
    % No FEC
    obj.setScenario('nofec');
    BER(2,i) = obj.simulate(EbN0_dB(i));
    % Hard dec
    obj.setScenario('convh');
    BER(3,i) = obj.simulate(EbN0_dB(i));
    % Soft dec
    obj.setScenario('convs');
    BER(4,i) = obj.simulate(EbN0_dB(i));
    % Puncturing!
    obj.setScenario('convp');
    BER(5,i) = obj.simulate(EbN0_dB(i));
    % Puncturing and RS!
    obj.setScenario('rsit0');
    BER(6,i) = obj.simulate(EbN0_dB(i));
    % Interleaving
    obj.setScenario('rsit1');
    BER(7,i) = obj.simulate(EbN0_dB(i));
    
    % Update Progress
    progressbar(i/length(EbN0_dB));
    
end
toc

% profile viewer

% Save values

save(['BER_Curves_' num2str(N_SYMBOLS)], 'EbN0_dB', 'BER');

% Plot BER Curves
h_ber_plot = figure;
semilogy(EbN0_dB,BER(1,:))
hold all
semilogy(EbN0_dB,BER(2,:))
semilogy(EbN0_dB,BER(3,:))
semilogy(EbN0_dB,BER(4,:))
semilogy(EbN0_dB,BER(5,:))
semilogy(EbN0_dB,BER(6,:))
semilogy(EbN0_dB,BER(7,:))

grid on
xlabel('E_b/N_0 [dB]')
ylabel('BER')
title('BER for Gray coded QPSK over AWGN')
legend('Theoretical', LABEL_SCENARIOS)

% Export figures...
saveas(h_ber_plot, ['BER_Curves_' num2str(N_SYMBOLS)], 'epsc')
saveas(h_ber_plot, ['BER_Curves_' num2str(N_SYMBOLS)], 'png')
