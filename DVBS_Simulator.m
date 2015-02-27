classdef DVBS_Simulator < handle
    %   DVBS_SIMULATOR Simulator class for BER calculation
    %   SCS Master SIMCOM DVB-S Simulator
    %   SCS 2014/2015 Juan Pablo Cuadro and Loic Veillard
    %
    % To set FEC options, class constructor can be used:
    %
    %       DVBS_Simulator(scenario_string, number_of_symbols)
    %
    % SCENARIO_STRING defines FEC options and it can also be 
    % set by using the setScenario method any time using the same
    % initialization string. Possible scenario strings:
    % 'nofec' : No FEC is used
    % 'convh' : Convolutional - 1/2 with hard decoding.
    % 'convs' : Convolutional - 1/2 with soft decoding.
    % 'convp' : Convolutional - 2/3 with puncturing.
    % 'rsit0' : RS + Convolutional - x/x with no interleaving.
    % 'rsit1' : RS + Convolutional - x/x with interleaving.
    %
    % NUMBER_OF_SYMBOLS defines the ammount of QPSK symbols transmitted in
    % the simulation. Bit padding in RS encoder is not implemented so care must be
    % taken in order to input a multiple of the RS encoder input bit count.
    
    properties  (Access = 'public')
        oversampling
        symbol_rate
        tx_filter
        coding
        dll
    end
    
    properties  (Access = 'public')
        simulation
    end
    
    methods     (Access = 'private')
        
        function this = resetScenario(this)
            % Convolutional coding on by default...
            this.coding.conv.switch = false;
            this.coding.conv.puncturingFlag = false;
            this.coding.conv.trellis = poly2trellis(7, [171, 133]);
            this.coding.conv.rate = 0.5;
            this.coding.conv.decType = 'hard';
            % RS coding off by default...
            this.coding.rs.switch = false;
            this.coding.rs.rate = 188/204;
            % Interleaving
            this.coding.interleaver.switch = false;
            this.coding.interleaver.depth = 10;
            this.dll.flag  = false;
            this.dll.order = 1;
            this.dll.d_phi_deg = 0;
            this.dll.BlT_dB = -3;
        end
                
        function this = tx(this)
            % FEC
            bits_encod = this.simulation.bit_stream_in;
            if this.coding.rs.switch %Reed-Solomon
                bits_encod = this.rs_enc(bits_encod);
            end
            if this.coding.interleaver.switch
                bits_encod = this.interleaver(bits_encod, this.coding.interleaver.depth);
            end            
            if this.coding.conv.switch %Convolutional
                bits_encod = this.conv_enc( bits_encod, this.coding.conv.trellis );
                % Puncturing
                if this.coding.conv.puncturingFlag
                    this.coding.conv.decType = 'soft';
                    [bits_encod, nPuncPad ] = this.puncture(bits_encod);
                    assert(nPuncPad == 0, 'Padded needed in puncture (not implemented)');
                end
            end
            % Pulse Shaping Filter
            this.tx_filter.h     = rcosfir(this.tx_filter.rolloff,...
                this.tx_filter.delay,...
                this.oversampling,...
                this.symbol_rate * this.oversampling / this.simulation.code_rate,...
                'sqrt');
            filt_delay  = this.tx_filter.delay;
            filt_h      = this.tx_filter.h;
            % Zero-padding
            bit_pad = zeros(1,2*2*filt_delay);
            bits_tx = [bits_encod, bit_pad]; % ZeroPadding
            % Modulation
            qpsk_symbols = this.qpskMapper(bits_tx);
            % Upsampling
            pad_v    = zeros(1,this.oversampling);
            pad_v(1) = 1;
            qpsk_symbols_over  = kron(qpsk_symbols, pad_v);
            % Filtering
            signal_tx = filter(filt_h, 1, qpsk_symbols_over);
            % Add signals to sink
            this.simulation.sink.bits_in  = this.simulation.bit_stream_in;
            this.simulation.sink.symbol_tx = qpsk_symbols;
            this.simulation.sink.symbol_tx_over = qpsk_symbols_over;
            this.simulation.sink.signal_tx = signal_tx;
            % DLL Init
            if this.dll.flag
                this.dll.pente = this.dll_pre(this.simulation.sink.symbol_tx_over,filt_h,this.oversampling);
            end
        end
        
        function this = ch(this,ebn0_db)
            signal_in = this.simulation.sink.signal_tx;
            EbN0 = 10^(ebn0_db/10);
            n_q_var = sum(abs(this.tx_filter.h).^2) * 1 / (4 * EbN0 / this.simulation.code_rate);
            noise = sqrt(n_q_var) * ( randn(1,length(signal_in)) + 1i*randn(1,length(signal_in)) );
            this.simulation.sink.signal_rx = signal_in + noise;
        end
        
        function this = rx(this)           
            % Matched filter
            signal_filt = filter(conj(fliplr(this.tx_filter.h)), 1, this.simulation.sink.signal_rx);
            % Sampling
            qpsk_symbols_hat = signal_filt(1:this.oversampling:end);
            % DLL
            if this.dll.flag
                [qpsk_symbols_hat, this.dll.debug] = this.dll_lock(...
                    qpsk_symbols_hat,...
                    this.dll.pente,...
                    this.dll.BlT_dB);
            end
            
            % Delay Compensation
            qpsk_symbols_hat = qpsk_symbols_hat((1+2*this.tx_filter.delay):end);
            bits_soft = reshape([real(qpsk_symbols_hat); imag(qpsk_symbols_hat)], 1, 2 * length(qpsk_symbols_hat));
            % FEC
            bits_decod = bits_soft < 0;
            if this.coding.conv.switch
                % Decode
                switch this.coding.conv.decType
                    case 'hard'
                        bits_decod = vitdec(bits_decod, this.coding.conv.trellis, 4, 'trunc', 'hard');
                    case 'soft'
                        if this.coding.conv.puncturingFlag
                            bits_soft = this.depuncture(bits_soft,0);
                        end
                        bits_decod = vitdec(bits_soft, this.coding.conv.trellis, 4, 'trunc', 'unquant');
                end
            end
            if this.coding.interleaver.switch
                bits_decod = this.deinterleaver(bits_decod, this.coding.interleaver.depth);
            end
            if this.coding.rs.switch
                bits_decod = this.rs_dec(bits_decod);
            end
            % Add signals to sink
            this.simulation.sink.bits_out = bits_decod;
            this.simulation.sink.symbol_rx = qpsk_symbols_hat;
            this.simulation.sink.signal_rx_amf = signal_filt;
        end
        
    end
    
    methods(Static)
        
        function pente = dll_pre(input,h,oversampling)
            d_phi_deg=[-3:3];
            for jj=1:length(d_phi_deg)  % boucle sur erreur de phase
                d_phi=d_phi_deg(jj)*pi/180+2*pi;
                input_dp=input.*exp(1i*d_phi);
                filtered = filter(h, 1, input_dp);
                decimate = filtered(1:oversampling:end);
                out_det = -imag(decimate.^4);
                S_curve(jj)=mean(out_det); 
            end
            pente=S_curve((length(S_curve)+1)/2+3)-S_curve((length(S_curve)+1)/2-3);
            pente=pente/(6*(d_phi_deg(2)-d_phi_deg(1))*pi/180);            
        end
        
        function [out, debug] = dll_lock(qpsk_symbols_hat, pente, BlT_dB)
                % Initialize DLL registers
                BlTT=10.^(BlT_dB);
                order = 1;
                if order==2
                    zeta=sqrt(2)/2;
                    A=16*zeta^2*BlTT.*(1+4*zeta^2-4*BlTT)/(1+4*zeta^2)./(1+4*zeta^2-8*zeta^2*BlTT);
                    B=64*zeta^2*BlTT.^2/(1+4*zeta^2)./(1+4*zeta^2-8*zeta^2*BlTT);
                elseif order==1
                    B=0*BlTT;
                    A=4*BlTT;
                end
                
                A=A/pente;
                B=B/pente;
                NCO_mem=0;      % initialisation NCO
                filtre_mem=0;   % initialisation de la memoire du filtre
                phi_est(1)=0;  % phase estimee : valeur initiale
                % DLL
                for ii=1:length(qpsk_symbols_hat)
                    out_det(ii)= -imag( (qpsk_symbols_hat(ii) * exp(-1i*phi_est(ii) )).^4 );
                    % filtre de boucle
                    w(ii)=filtre_mem+out_det(ii); % memoire filtre + sortie detecteur
                    filtre_mem=w(ii);
                    out_filtre=A*out_det(ii)+B*w(ii);   % sortie du filtre a l'instant ii :  F(z)=A+B/(1-z^-1)
                    %NCO
                    phi_est(ii+1)=(out_filtre+NCO_mem); % N(z)=1/(z-1)
                    NCO_mem=phi_est(ii+1);
                end
                out = qpsk_symbols_hat .* exp(-1i*phi_est(2:end));
                
                debug.phi_est = phi_est;
        end
        
        function [b_out] = conv_enc( b_in, trellis )
            b_out = convenc( b_in, trellis );
        end
        
        function [b_out, P] = puncture(b_in)
            punc_vec = logical([1 1 0 1]);
            L = length(punc_vec);
            N = length(b_in);
            P = mod(L - mod(N, L), L);
            b_out = [b_in zeros(1, P)];
            mask = repmat(~punc_vec, 1, (N+P)/L);
            b_out(mask) = [];
        end
        
        function [b_soft_out] = depuncture(b_soft_in, P)
            punc_vec = logical([1 1 0 1]);
            M = length(b_soft_in) / sum(punc_vec);
            mask  = repmat(punc_vec, 1, M);
            b_soft_out = zeros(1, M * length(punc_vec));
            b_soft_out(mask) = b_soft_in;
            b_soft_out = b_soft_out(1:(end-P));
        end
        
        function output = interleaver(input, block_interleave_depth)
            nbits = 8;
            rs_codeword_size_symbols = 204;

            temp_binary = reshape(input, nbits, []).';            
            temp_decima = temp_binary * (2.^(size(temp_binary,2) - 1:-1:0))';
            temp_symbol = reshape(temp_decima,...
                rs_codeword_size_symbols * block_interleave_depth, []);
            
            output = nan(size(temp_symbol));
            for i = 1:size(output,2)
                output(:,i) = matintrlv(temp_symbol(:,i), ...
                    block_interleave_depth, rs_codeword_size_symbols);
            end
            
            output = output(:);
            output = (dec2bin(output, nbits)-'0').';
            output = output(:).';
        end
        
        function output = deinterleaver(input, block_interleave_depth)
            nbits = 8;
            rs_codeword_size_symbols = 204;

            temp_binary = reshape(input, nbits, []).';            
            temp_decima = temp_binary * (2.^(size(temp_binary,2) - 1:-1:0))';
            temp_symbol = reshape(temp_decima,...
                rs_codeword_size_symbols * block_interleave_depth, []);
            
            output = nan(size(temp_symbol));
            for i = 1:size(output,2)
                output(:,i) = matdeintrlv(temp_symbol(:,i), ...
                    block_interleave_depth, rs_codeword_size_symbols);
            end
            output = output(:);
            output = (dec2bin(output, nbits)-'0').';
            output = output(:).';
        end
        
        function c_out = qpskMapper(b_in)
            Nb = length(b_in);
            Ns = Nb/2;
            IQ = (1-reshape(b_in,[2,Ns]).*2) * cos(pi/4);
            c_out = IQ(1,:) + 1i*IQ(2,:);
        end
        
        function b_out = rs_enc(b_in)
            warning('off')
            M = 8;
            S = 51; % Shortening
            t = 8;
            Ks = 188;
            assert(mod(length(b_in), Ks*M) == 0,'bit_stream length must be a multiple of 1504');
            % Add 51 bytes (set to 0)
            padded = reshape(b_in.', Ks*M, []).';
            padded = [zeros(size(padded,1), M*S), padded];
            padded = reshape(padded.', length(padded(:)), []).';
            %RS(255,203,8) encoding           
            N=2^M-1;   % nombre de symboles du mot de code RS (apres codage)
            K=N-2*t;   % nombre de symboles du mot d'info RS
            nb_mots=length(padded)/(M*K);
            %
            % table gf2ded pour transformer un symbole non nul m d'un corps de Gallois en un
            % nombre non nul m (attention, m=0 doit etre traiter a part)
            %            
            gf_nonzero=gf([1:2^M-1],M);
            expformat=log(gf_nonzero)+1;
            for ii=1:2^M-1
                table_gf2dec(expformat(ii))=ii;
            end
            %
            % codage
            %
            mess=padded;   % bits correspondants a (nb_mots) mots de K*M bits (vecteur)
            mess_matrix=reshape(mess,M,K*nb_mots)'; % matrice de (nb_mots*K) colonnes de M bits
            mess_matrix_symb=bi2de(mess_matrix);       % matrice de nb_mots mots de K entiers b (0<=b<=2^M-1)
            mots_info=reshape(mess_matrix_symb',1,nb_mots*K); %vecteurs de nb_mots d'info de K entiers b
            msg  = gf(mots_info,M);                           % les entiers b sont transformes en elements de GF
            msg_forme=reshape(msg,K,nb_mots)';
            mots_code=rsenc(msg_forme,N,K);   %#ok<REMFF1> %codage : on obtient N elements de GF en sortie
            mots_code_serie=reshape(mots_code',1,N*nb_mots);
            mots_code_de=DVBS_Simulator.gf2dec(mots_code_serie,table_gf2dec);  %  symbole de Gf => entiers
            mat_bits=de2bi(mots_code_de');                % passage entiers => bits (matrice)
            b_out=reshape(mat_bits',1,nb_mots*N*M);   % pasage matrice bits => vecteurs de bits
            
            % Shortening
            b_out = reshape(b_out.', N*M, []).';
            b_out(:,1:(S*M)) = [];
            b_out = reshape(b_out.', length(b_out(:)), []).';
            warning('on')
        end
        
        function b_out = rs_dec(b_in)
            warning('off')
            M = 8;
            S = 51; % Shortening
            N = 2^M-1;
            Ks = 188;
            K  = Ks + S;
            % Add 51 bytes (set to 0)
            padded = reshape(b_in.', (N-S)*M, []).';
            padded = [zeros(size(padded,1), M*S), padded];
            padded = reshape(padded.', length(padded(:)), []).';
            %RS(255,203,8) decoding
            nb_mots=length(padded)/(N*M);
            %
            % table gf2ded pour transformer un symbole non nul m d'un corps de Gallois en un
            % nombre non nul m (attention, m=0 doit etre traiter a part)
            %
            gf_nonzero=gf([1:2^M-1],M);
            expformat=log(gf_nonzero)+1;
            for ii=1:2^M-1
                table_gf2dec(expformat(ii))=ii;
            end
            bits_code_matrix=reshape(padded,M,N*nb_mots)'; % matrice de (nb_mots*N) colonnes de M bits
            symb_matrix=bi2de(bits_code_matrix);       % matrice de nb_mots mots de N entiers b (0<=b<=2^M-1)
            mots_recu=reshape(symb_matrix',1,nb_mots*N); %vecteurs de nb_mots d'info de N entiers b
            msg_recu  = gf(mots_recu,M);  % les entiers b sont transformes en elements de GF
            msg_recu_forme=reshape(msg_recu,N,nb_mots)';
            mots_decode=rsdec(msg_recu_forme,N,K);   %decodage : on obtient K elements de GF en sortie
            mots_decode_serie=reshape(mots_decode',1,K*nb_mots);
            mots_decode_de=DVBS_Simulator.gf2dec(mots_decode_serie,table_gf2dec);  %  symbole de Gf => entiers
            mat_bits_recu=de2bi(mots_decode_de');                % passage entiers => bits (matrice)
            b_out=reshape(mat_bits_recu',1,nb_mots*K*M);   % pasage matrice bits => vecteurs de bits

            % Shortening
            b_out = reshape(b_out.', K*M, []).';
            b_out(:,1:(S*M)) = [];
            b_out = reshape(b_out.', length(b_out(:)), []).';
            warning('on')
        end
        
        function a_dec=gf2dec(a,tabb)
            a_dec(1:length(a))=0;
            nzg = find(a~=0);
            a_dec(nzg)=tabb(log(a(nzg))+1);
        end
        
    end
    
    methods     (Access = 'public')
        
        function this = DVBS_Simulator(scenario_string, number_of_symbols)
            % Default values...
            this.oversampling = 10;
            this.symbol_rate = 1e3;
            this.tx_filter.delay = 6;
            this.tx_filter.rolloff = 0.35;
            
            % FEC Options
            this.setScenario(scenario_string);
            
            % Initialize bit_stream
            n_bits = 2 * number_of_symbols;
            this.simulation.bit_stream_in = round(rand(1,n_bits));
        end
        
        function this = setScenario(this, string)
            % Possible scenarios:
            % 'nofec' : No FEC is used
            % 'convh' : Convolutional - 1/2 with hard decoding.
            % 'convs' : Convolutional - 1/2 with soft decoding.
            % 'convp' : Convolutional - 2/3 with puncturing.
            % 'rsit0' : Reed Solomon enabled and interleaving on.
            % 'rsit1' : Reed Solomon enabled and interleaving off.
            % 'rsonl' : Reed Solomon enabled and convolutional disabled.
            
            this.resetScenario();            
            
            switch string
                case 'nofec'
                    this.coding.conv.switch = false;
                    this.coding.conv.puncturingFlag = false;
                    this.coding.conv.decType = 'hard';
                    this.coding.rs.switch = false;
                    this.coding.interleaver.switch = false;
                    this.coding.interleaver.depth = 10;
                case 'convh'
                    this.coding.conv.switch = true;
                    this.coding.conv.puncturingFlag = false;
                    this.coding.conv.decType = 'hard';
                    this.coding.rs.switch = false;
                    this.coding.interleaver.switch = false;
                    this.coding.interleaver.depth = 10;
                case 'convs'
                    this.coding.conv.switch = true;
                    this.coding.conv.puncturingFlag = false;
                    this.coding.conv.decType = 'soft';
                    this.coding.rs.switch = false;
                    this.coding.interleaver.switch = false;
                    this.coding.interleaver.depth = 10;
                case 'convp'
                    this.coding.conv.switch = true;
                    this.coding.conv.puncturingFlag = true;
                    this.coding.conv.decType = 'soft';
                    this.coding.rs.switch = false;
                    this.coding.interleaver.switch = false;
                    this.coding.interleaver.depth = 10;
                case 'rsit0'
                    this.coding.conv.switch = true;
                    this.coding.conv.puncturingFlag = true;
                    this.coding.conv.decType = 'soft';
                    this.coding.rs.switch = true;
                    this.coding.interleaver.switch = false;
                    this.coding.interleaver.depth = 10;
                case 'rsit1'
                    this.coding.conv.switch = true;
                    this.coding.conv.puncturingFlag = true;
                    this.coding.conv.decType = 'soft';
                    this.coding.rs.switch = true;
                    this.coding.interleaver.switch = true;
                    this.coding.interleaver.depth = 10;
                case 'rsonl'
                    this.coding.conv.switch = false;
                    this.coding.conv.puncturingFlag = false;
                    this.coding.conv.decType = 'hard';
                    this.coding.rs.switch = true;
                    this.coding.interleaver.switch = false;
                    this.coding.interleaver.depth = 10;
                otherwise
                    error('Scenario not valid check help.');
            end
                                
        end
        
        function ber = simulate(this,ebn0_db)
            % Calculate Code Rate
            code_rate = 1;
            if this.coding.conv.switch
                code_rate = this.coding.conv.rate;
                if this.coding.conv.puncturingFlag
                    code_rate = 2/3;                     
                end
            end
            if this.coding.rs.switch
                code_rate = code_rate * this.coding.rs.rate;
            end
            this.simulation.code_rate = code_rate;
            % Tx and Rx
            this = this.tx();
            this = this.ch(ebn0_db);
            this = this.rx();
            % Calculate BER
            this.simulation.ber = sum(this.simulation.bit_stream_in ~=...
                this.simulation.sink.bits_out) / length(this.simulation.sink.bits_out);
            ber = this.simulation.ber;
        end
        
        function figure_handle = plotConstellation(this)
            ang=0:0.01:2*pi;
            xp=cos(ang);
            yp=sin(ang);
                    
            
            % Plot
            figure_handle = figure;
            hold all
            
            h1 = plot(this.simulation.sink.symbol_rx,'.','MarkerSize',5);
            plot(xp,yp,':')
            h2 = plot(this.simulation.sink.symbol_tx,'or','MarkerSize', 6);
            set(h1, 'MarkerFaceColor', get(h1, 'Color'));
            set(h2, 'MarkerFaceColor', get(h2, 'Color'));
            xlim([-1.5 1.5])
            ylim([-1.5 1.5])
            xlabel('In-phase Amplitude')
            ylabel('Quadrature Amplitude')
            title('IQ samples')
            box on
            legend([h1,h2],'Rx','Tx','Location','EastOutside')
            line([-1.5 1.5],[0 0],'Color','k')
            line([0 0],[-1.5 1.5],'Color','k')
            grid on
            
        end
        
        function figure_handle = plotSignals(this)
            Ts   = 1/(this.symbol_rate / this.simulation.code_rate);
            N_up = this.oversampling;
            signal_tx = this.simulation.sink.signal_tx;
            signal_rx = this.simulation.sink.signal_rx_amf;
            timebase_s = (0:length(signal_tx)-1) * Ts / N_up ;
            timebase_Ts = timebase_s / Ts;
            % Plot
            figure_handle = figure;
            h1= subplot(2,1,1);
            hold all
            plot(timebase_Ts, real(signal_tx), 'LineWidth', 2)
            plot(timebase_Ts, real(signal_rx))
            xlabel('Time/Ts')
            title('In-Phase Signal Components')
            h2 = subplot(2,1,2);
            hold all
            plot(timebase_Ts, imag(signal_tx), 'LineWidth', 2)
            plot(timebase_Ts, imag(signal_rx))
            xlabel('Time/Ts')
            title('Quadrature Signal Components')
            linkaxes([h1,h2],'xy');
            xlim([0,timebase_Ts(end)])            
        end
        
        function figure_handle = plotPsd(this)
            Rs   = this.symbol_rate / this.simulation.code_rate;
            N_up = this.oversampling;
            signal_tx = this.simulation.sink.signal_tx;
            signal_rx = this.simulation.sink.signal_rx;
            % Plot
            figure_handle = figure;
            [psd_tx,F] = pwelch(signal_tx,[],[],[],Rs*N_up,'centered');
            [psd_rx,~] = pwelch(signal_rx,[],[],[],Rs*N_up,'centered');
            F_norm = F / Rs;
            hold on
            plot(F_norm, 10*log10(psd_tx))
%             plot(F_norm, 10*log10(psd_rx),'r')
            title('Power Spectral Density')
            xlabel('Frequency')
            ylabel('Power Spectral Density [dBW/Hz]')
%             legend('Tx','Rx')
            grid on
            box on
            xlim([-1,1]);
            set(gca,'XTick',[-1, -0.5, 0, 0.5, 1])
            set(gca,'XTickLabel',{'-R_s', '-R_s/2', '0', 'R_s/2', 'R_s'})
        end
        
        function figure_handle = plotEyes(this, periods)
            N_up = this.oversampling;
            signal_tx = this.simulation.sink.signal_tx;
            signal_rx = this.simulation.sink.signal_rx_amf;
            
            idx = 1:(periods * N_up);
            % Plot
            figure_handle = eyediagram(signal_rx(20*N_up + idx), N_up, 1, 0, 'k');
            figure_handle.Children(1).XLabel = xlabel('Time / Ts');
            figure_handle.Children(2).XLabel = xlabel('Time / Ts');
        end
        
    end
    
end

