pente1 = PLL_QPSK_BO_NDA([-180:180],2);
export_fig 2-1_1_007dB.pdf -transparent
close
pente2 = PLL_QPSK_BO_NDA([-180:180],1e16);
export_fig 2-1_1_---dB.pdf -transparent
close


%% Plot against BLT

d_phi_deg = 10;
df_Rs = 0;
order = 2;
BlT_dB = [-1, -2, -3, -4];
EbNodB = 100;
clear phi_est_deg
clear B_w
clear out_det

for ii = 1:numel(BlT_dB)
   [phi_est_deg(ii,:), B_w(ii,:), out_det(ii,:)] = Pll_qpsk_NDA(d_phi_deg, df_Rs, BlT_dB(ii), order, EbNodB);
end

figure
hold on
plot(phi_est_deg.')
plot([1,size(phi_est_deg,2)], [d_phi_deg,d_phi_deg], '-k')
grid on
xlabel('Time')
ylabel('Phi Estimate [degrees]')
legend('B_LT = 10^{-1}','B_LT = 10^{-2}', 'B_LT = 10^{-3}', 'B_LT = 10^{-4}');
title('E_b/N_0 = 100dB') 

export_fig 2-2_1_100dB.pdf -transparent
close

%% Plot against BLT

d_phi_deg = 10;
df_Rs = 0;
order = 1;
BlT_dB = [-1, -2, -3, -4];
EbNodB = 7;
clear phi_est_deg
clear B_w
clear out_det

for ii = 1:numel(BlT_dB)
   [phi_est_deg(ii,:), B_w(ii,:), out_det(ii,:)] = Pll_qpsk_NDA(d_phi_deg, df_Rs, BlT_dB(ii), order, EbNodB);
end

figure
hold on
plot(phi_est_deg.')
plot([1,size(phi_est_deg,2)], [d_phi_deg,d_phi_deg], '-k')
grid on
xlabel('Time')
ylabel('Phi Estimate [degrees]')
legend('B_LT = 10^{-1}','B_LT = 10^{-2}', 'B_LT = 10^{-3}', 'B_LT = 10^{-4}');
title('E_b/N_0 = 7dB') 

export_fig 2-2_1_007dB.pdf -transparent
close

%% Plot against BLT

d_phi_deg = 30;
df_Rs = 0.01;
order = 2;
BlT_dB = [-1, -2, -3];
EbNodB = 100;
clear phi_est_deg
clear B_w
clear out_det

for ii = 1:numel(BlT_dB)
   [phi_est_deg(ii,:), B_w(ii,:), out_det(ii,:)] = Pll_qpsk_NDA(d_phi_deg, df_Rs, BlT_dB(ii), order, EbNodB);
end

figure
plot(B_w.')
grid on
xlabel('Time')
ylabel('Frequency error');
legend('B_LT = 10^{-1}','B_LT = 10^{-2}', 'B_LT = 10^{-3}');
title('E_b/N_0 = 100dB') 

export_fig 2-2_2_100dB.pdf -transparent
close
%% Plot against BLT

d_phi_deg = 30;
df_Rs = 0.01;
order = 2;
BlT_dB = [-1, -2, -3];
EbNodB = 7;
clear phi_est_deg
clear B_w
clear out_det

for ii = 1:numel(BlT_dB)
   [phi_est_deg(ii,:), B_w(ii,:), out_det(ii,:)] = Pll_qpsk_NDA(d_phi_deg, df_Rs, BlT_dB(ii), order, EbNodB);
end

figure
plot(B_w.')
grid on
xlabel('Time')
ylabel('Frequency error');
legend('B_LT = 10^{-1}','B_LT = 10^{-2}', 'B_LT = 10^{-3}');
title('E_b/N_0 = 7dB') 

export_fig 2-2_2_007dB.pdf -transparent
close