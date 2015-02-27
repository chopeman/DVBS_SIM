%   SIMCOM DVB-S Simulator
%   2014/2015 Juan Pablo Cuadro and Loic Veillard

close all
clear 

cd report
delete *.eps

cd ..

% PSD
N_SYMBOLS = 20000;
sim = DVBS_Simulator('nofec', N_SYMBOLS);
sim.simulate(10);
hpsd = sim.plotPsd;
set(hpsd,'PaperUnits','points','PaperPosition',[300 400 1.6*560 420])
saveas(hpsd, 'report/matp_TxPsd', 'epsc')
close(hpsd)

% Constellation
N_SYMBOLS = 1000;
EbN0_dB = 10;
sim = DVBS_Simulator('nofec', N_SYMBOLS);
sim.simulate(EbN0_dB);
hpsd = sim.plotConstellation;
print('report/matp_Constellation10dB', '-dpng','-r300')
close(hpsd)

% Eye Diagrams
N_SYMBOLS = 1000;
sim = DVBS_Simulator('nofec', N_SYMBOLS);
sim.simulate(90);
hpsd = sim.plotEyes(200);
saveas(hpsd, 'report/matp_EyeClean', 'epsc')
close(hpsd)

sim.simulate(10);
hpsd = sim.plotEyes(200);
saveas(hpsd, 'report/matp_Eye10dB', 'epsc')
close(hpsd)

% Load BER values
load BER_Curves_300800.mat;

labels = {'No FEC, Theoretical',...
    'No FEC, Sim',...
    'Convolutional Hard Decoder 1/2, Sim',...
    'Convolutional Soft Decoder 1/2, Sim',...
    'Convolutional Punctured 2/3, Sim',...
    'Reed Solomon 2/3, Sim',...
    'Reed Solomon Interleaved 2/3, Sim',...  
    'Reed Solomon Only, Sim'
    };
% Plot BER Curves
c = [0, 113, 188; 125, 46, 141; 216, 82, 24;]/255;
f = figure;
semilogy(EbN0_dB,BER(1,:), 'color', c(1,:))
hold all
semilogy(EbN0_dB,BER(2,:),'diamond', 'color', c(1,:))

semilogy(EbN0_dB,BER(3,:), '-', 'color', c(2,:))
semilogy(EbN0_dB,BER(4,:), '--', 'color', c(2,:))
semilogy(EbN0_dB,BER(5,:), '-.', 'color', c(2,:))

semilogy(EbN0_dB,BER(6,:), 'color', c(3,:))
semilogy(EbN0_dB,BER(7,:), '--', 'color', c(3,:))
semilogy(EbN0_dB,BER(8,:), '-.', 'color', c(3,:))
grid on
xlabel('E_b/N_0 [dB]')
ylabel('BER')
title('BER for Gray coded QPSK over AWGN')
legend(labels)

set(f, 'units', 'inches', 'position', [ 0 0 10 7])
set(f,'PaperUnits','inches','PaperPosition',[0 0 12 9])
saveas(f, 'report/matp_BER', 'epsc')
close(f)




