clear
clc
close all

N_SYMBOLS = 1e5;
SCENARIOS = [-2, -2.5, -3, -3.5, -4];
LABEL_SCENARIOS = {'Theoretical',...
    'log(BT) = -2',...
    'log(BT) = -2.5',...
    'log(BT) = -3',...
    'log(BT) = -3.5',...
    'log(BT) = -4';
    };
N_SCENARIOS  = numel(LABEL_SCENARIOS);

EbN0_dB = -5:3:10;
EbN0_dB_2 = EbN0_dB(1):1:EbN0_dB(end);
BER_2 = 0.5*erfc(sqrt(10.^(EbN0_dB_2/10)));

BER = zeros(N_SCENARIOS, length(EbN0_dB));

% Theoretical BER
BER(1,:) = 0.5*erfc(sqrt(10.^(EbN0_dB/10)));

fprintf('\nStarting simulation...\n');
fprintf('Number of bits = %d\n',2*N_SYMBOLS);
fprintf('Number scenarios = %d\n',N_SCENARIOS);
fprintf('Number EbN0 points = %d\n\n',length(EbN0_dB));
fprintf('\n\n');
tic

loop_counter = 0;
loop_total = 3 * length(EbN0_dB);


% Construct Simulator Class
obj = DVBS_Simulator('nofec', N_SYMBOLS);

phi_est_cell = cell(length(SCENARIOS),length(EbN0_dB));

for ii = 1:length(SCENARIOS)
    obj.dll.BlT_dB = SCENARIOS(ii);
    obj.dll.d_phi_deg = 0;
    obj.dll.flag = true;
    for jj = 1:length(EbN0_dB)
        BER(ii+1,jj) = obj.simulate(EbN0_dB(jj));
        phi_est_cell{ii,jj} = obj.dll.debug.phi_est;
        loop_counter = loop_counter + 1;
        progressDisplay(loop_counter,loop_total);
    end
    
end

close all
set(0, 'defaultTextInterpreter', 'latex');

% Plot BER Curves
f = figure;
semilogy(EbN0_dB_2,BER_2)
hold all
semilogy(EbN0_dB,BER(2,:))
semilogy(EbN0_dB,BER(3,:))
semilogy(EbN0_dB,BER(4,:))
semilogy(EbN0_dB,BER(5,:))
semilogy(EbN0_dB,BER(6,:))

grid on
xlabel('$$E_b/N_0 (dB)$$')
ylabel('BER')
legend(LABEL_SCENARIOS)
legend('Location','SouthWest');
export_fig DLL1.png -r300

figure
hold all
plot(rad2deg(phi_est_cell{2,5}))
plot(rad2deg(phi_est_cell{3,5}))
plot(rad2deg(phi_est_cell{4,5}))
grid on
ylim([-200,200])
xlim([0,1e5])
legend(LABEL_SCENARIOS(2:4))
xlabel('Symbol Index')
ylabel('$$\tilde{\theta} (^{\circ})$$')
export_fig DLL2.png -r300






