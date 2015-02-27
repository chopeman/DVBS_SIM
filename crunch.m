function crunch(order)
%   SIMCOM DVB-S Simulator
%   2014/2015 Juan Pablo Cuadro and Loic Veillard
% add progress display to end as subfunction!

N_SYMBOLS = 4 * 5 * 1.504 * order;
LABEL_SCENARIOS = {'Theoretical',...
    'Simulated - no FEC',...
    'Simulated - Conv Hard Deco - R = 1/2',...
    'Simulated - Conv Soft Deco - R = 1/2',...
    'Simulated - Conv Punctured - R = 2/3',...
    'Simulated - No Interleave  - RS',...
    'Simulated - Interleaved    - RS',...
    'Simulated - RS only'
    };
N_SCENARIOS     = numel(LABEL_SCENARIOS);

EbN0_dB = -4:1:7;

BER = zeros(N_SCENARIOS, length(EbN0_dB));

% Construct Simulator Class
obj = DVBS_Simulator('nofec', N_SYMBOLS);

% Theoretical BER
BER(1,:) = 0.5*erfc(sqrt(10.^(EbN0_dB/10)));

fprintf('\nStarting simulation...\n');
fprintf('Number of bits = %d\n',2*N_SYMBOLS);
fprintf('Number scenarios = %d\n',N_SCENARIOS);
fprintf('Number EbN0 points = %d\n\n',length(EbN0_dB));

tic

loop_counter = 0;
loop_total = 7 * length(EbN0_dB);

for i = 1:length(EbN0_dB)
    
    % **** Scenarios! **** %
    % No FEC
    obj.setScenario('nofec');
    BER(2,i) = obj.simulate(EbN0_dB(i));
    loop_counter = loop_counter + 1;
    progressDisplay(loop_counter,loop_total);
    % Hard dec
    obj.setScenario('convh');
    BER(3,i) = obj.simulate(EbN0_dB(i));
    loop_counter = loop_counter + 1;
    progressDisplay(loop_counter,loop_total);
    % Soft dec
    obj.setScenario('convs');
    BER(4,i) = obj.simulate(EbN0_dB(i));
    loop_counter = loop_counter + 1;
    progressDisplay(loop_counter,loop_total);
    % Puncturing!
    obj.setScenario('convp');
    BER(5,i) = obj.simulate(EbN0_dB(i));
    loop_counter = loop_counter + 1;
    progressDisplay(loop_counter,loop_total);
    
    % Puncturing and RS!
    obj.setScenario('rsit0');
    BER(6,i) = obj.simulate(EbN0_dB(i));
    loop_counter = loop_counter + 1;
    progressDisplay(loop_counter,loop_total);
    % Interleaving
    obj.setScenario('rsit1');
    BER(7,i) = obj.simulate(EbN0_dB(i));
    loop_counter = loop_counter + 1;
    progressDisplay(loop_counter,loop_total);
    % Only RS
    obj.setScenario('rsonl');
    BER(8,i) = obj.simulate(EbN0_dB(i));
    loop_counter = loop_counter + 1;
    progressDisplay(loop_counter, loop_total);
    
end

fprintf('\n\n');
toc

% Save values

save(['BER_Curves_' num2str(N_SYMBOLS)], 'EbN0_dB', 'BER');

end

function progressDisplay( iter, total )
if iter > 1;
    fprintf('\b\b\b')
elseif iter == total
    return
end
fprintf('%2d%s', floor(100*iter/total), '%');
end
