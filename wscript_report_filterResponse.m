close all
clear 

sim = DVBS_Simulator('nofec',1e3);
sim.simulate(0);
h = sim.tx_filter.h;
tap = 0.1 * (1:length(sim.tx_filter.h))-0.1;
freq = -2:0.001:2;
tr = nfourier2( tap, sim.tx_filter.h, freq)/sqrt(10);

fh = figure;
% plot(freq, 20*log10(abs(tr)))
plot(freq, (abs(tr).^2))
grid on
xlim([-1,1])
ylim([0,1.1])
% ylim([-60,1])
xlabel('Frequency / R_s')
ylabel('|H(f)|^2')
title('Frequency response of emitter filter');

set(gca,'YTick',0:0.1:1)
set(gca,'XTick',[-1, -1.35/2, -0.5, 0, 0.5, 1.35/2,  1])
set(fh,'PaperUnits','points','PaperPosition',[300 400 1.6*560 420])
saveas(fh, 'report/matp_txfiltresp', 'epsc')
close(fh)