%   2014/2015 Juan Pablo Cuadro

close all
clear

N_SYMBOLS = 4 * 5 * 1.504 * 1e4;
OVERSAMPLING = 10;
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

bers = zeros(N_SCENARIOS, length(EbN0_dB));

% Construct Simulator Class
obj = DVBS_Simulator(N_SYMBOLS);
obj.oversampling = OVERSAMPLING;

tic
bers(1,:) = 0.5*erfc(sqrt(10.^(EbN0_dB/10)));
for i = 1:length(EbN0_dB)
    progressbar(i/length(EbN0_dB));
    % No FEC
    obj.coding.conv.switch          = false;
    obj.coding.rs.switch            = false;
    obj.coding.conv.puncturingFlag  = false;
    bers(2,i) = obj.simulate(EbN0_dB(i));
    % Hard dec
    obj.coding.conv.switch          = true;
    obj.coding.conv.decType         = 'hard';
    obj.coding.rs.switch            = false;
    obj.coding.conv.puncturingFlag  = false;
    bers(3,i) = obj.simulate(EbN0_dB(i));
    % Soft dec
    obj.coding.conv.switch          = true;
    obj.coding.conv.decType         = 'soft';
    obj.coding.rs.switch            = false;
    obj.coding.conv.puncturingFlag  = false;
    bers(4,i) = obj.simulate(EbN0_dB(i));
    % Puncturing!
    obj.coding.conv.switch          = true;
    obj.coding.conv.decType         = 'soft';
    obj.coding.rs.switch            = false;
    obj.coding.conv.puncturingFlag  = true;
    bers(5,i) = obj.simulate(EbN0_dB(i));
    % Puncturing and RS!
    obj.coding.conv.switch          = true;
    obj.coding.conv.decType         = 'hard';
    obj.coding.rs.switch            = true;
    obj.coding.conv.puncturingFlag  = true;
    bers(6,i) = obj.simulate(EbN0_dB(i));
    % Interleaving
    obj.coding.conv.switch          = true;
    obj.coding.conv.decType         = 'hard';
    obj.coding.rs.switch            = true;
    obj.coding.conv.puncturingFlag  = true;
    obj.coding.interleaving         = true;
    bers(7,i) = obj.simulate(EbN0_dB(i));
    obj.coding.interleaving         = false;
    
end
toc

% Plot BER Curves

h_ber_plot = figure;
semilogy(EbN0_dB,bers(1,:))
hold all
semilogy(EbN0_dB,bers(2,:))
semilogy(EbN0_dB,bers(3,:))
semilogy(EbN0_dB,bers(4,:))
semilogy(EbN0_dB,bers(5,:))
semilogy(EbN0_dB,bers(6,:))
semilogy(EbN0_dB,bers(7,:))

grid on
xlabel('E_b/N_0 [dB]')
ylabel('BER')
title('BER for Gray coded QPSK over AWGN')
legend('Theoretical', LABEL_SCENARIOS)

saveas(h_ber_plot, ['BER_Curves_' num2str(N_SYMBOLS)], 'epsc')
saveas(h_ber_plot, ['BER_Curves_' num2str(N_SYMBOLS)], 'png')
