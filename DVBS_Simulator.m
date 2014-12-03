
classdef DVBS_Simulator < handle
    %DVBS_SIMULATOR Simulator class for DVB-S Tx and Rx
    %   Detailed explanation goes here
    
    properties (Access = 'public')
        oversampling
        symbol_rate
        tx_filter
        coding
    end
    
    properties (Access = 'private')
        simulation
    end
    
    methods (Access = 'private')
        
        function this = tx(this)
            % FEC
            bits_encod = this.simulation.bit_stream_in;
            if this.coding.rs.switch %Reed-Solomon
                bits_encod = this.rs_enc(bits_encod);
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
            this.simulation.sink.bits_FECenc_in  = bits_tx;
            this.simulation.sink.bits_FECenc_out = bits_encod;
            this.simulation.sink.symbol_tx = qpsk_symbols;
            this.simulation.sink.signal_tx = signal_tx;
            
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
            if this.coding.rs.switch
                bits_decod = this.rs_dec(bits_decod);
            end
            % Add signals to sink
            this.simulation.sink.bits_FECdec_out = bits_decod;
            this.simulation.sink.symbol_rx = qpsk_symbols_hat;
            
        end
        
    end
    
    methods(Static)
        
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
    
    
    % Public methods
    methods (Access = 'public')
        function this = DVBS_Simulator(varargin)
            % Default values...
            this.oversampling = 10;
            this.symbol_rate = 1e3;
            this.tx_filter.delay = 6;
            this.tx_filter.rolloff = 0.35;
            % Convolutional coding on by default...
            this.coding.conv.switch = true;
            this.coding.conv.puncturingFlag = false;
            this.coding.conv.trellis = poly2trellis(7, [171, 133]);
            this.coding.conv.rate = 0.5;
            this.coding.conv.decType = 'hard';
            % RS coding off by default...
            this.coding.rs.switch = false;
            this.coding.rs.rate = 1;
            % Initialize bit_stream
            if nargin > 0
                n_bits = 2*round(varargin{1});
            else
                n_bits = 2*1e3;
            end
            this.simulation.bit_stream_in = round(rand(1,n_bits));
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
                this.simulation.sink.bits_FECdec_out) / length(this.simulation.sink.bits_FECdec_out);
            ber = this.simulation.ber;
        end
        
        function plotConstellation(this)
            % Plot
            figure
            hold all
            plot(this.simulation.sink.symbol_rx,'.b','MarkerSize',5)
            plot(this.simulation.sink.symbol_tx, 'or','MarkerSize', 6, 'MarkerFaceColor', 'r')
            xlim([-1.5 1.5])
            ylim([-1.5 1.5])
            xlabel('In-phase Amplitude')
            ylabel('Quadrature Amplitude')
            title('QPSK Constellation')
            box on
            legend('Rx','Tx','Location','EastOutside')
            line([-1.5 1.5],[0 0],'Color','k')
            line([0 0],[-1.5 1.5],'Color','k')
            grid on
            
        end
        
        function plotSignals(this)
            Rs   = this.symbol_rate / this.simulation.code_rate;
            N_up = this.oversampling;
            signal_tx = this.simulation.sink.signal_tx;
            signal_rx = this.simulation.sink.signal_rx;
            timebase = (0:length(signal_tx)-1) * Rs * N_up ;
            % Plot
            figure
            h1= subplot(3,1,1);
            hold all
            plot(timebase, real(signal_tx), 'b')
            plot(timebase, real(signal_rx), 'r')
            xlabel('Time [s]')
            title('In-Phase Signal Components')
            h2 = subplot(3,1,2);
            hold all
            plot(timebase, imag(signal_tx), 'b')
            plot(timebase, imag(signal_rx), 'r')
            xlabel('Time [s]')
            title('Quadrature Signal Components')
            linkaxes([h1,h2],'xy');
            xlim([0,timebase(end)])
            subplot(3,1,3)
            [psd_tx,F] = pwelch(signal_tx,[],[],[],Rs*N_up,'centered');
            [psd_rx,~] = pwelch(signal_rx,[],[],[],Rs*N_up,'centered');
            hold on
            plot(F, 10*log10(psd_tx),'b')
            plot(F, 10*log10(psd_rx),'r')
            title('Power Spectral Density')
            
        end
        
        function plotPsd(this)
            Rs   = this.symbol_rate / this.simulation.code_rate;
            N_up = this.oversampling;
            signal_tx = this.simulation.sink.signal_tx;
            signal_rx = this.simulation.sink.signal_rx;
            % Plot
            figure
            [psd_tx,F] = pwelch(signal_tx,[],[],[],Rs*N_up/1e3,'centered');
            [psd_rx,~] = pwelch(signal_rx,[],[],[],Rs*N_up/1e3,'centered');
            hold on
            plot(F, 10*log10(psd_tx),'b')
            plot(F, 10*log10(psd_rx),'r')
            title('Power Spectral Density')
            xlabel('Frequency [kHz]')
            ylabel('Power Spectral Density [W/Hz]')
            legend('Tx','Rx')
            grid on
            box on
        end
        
    end
    
end

