function [phi_est_deg, B_w, out_det] = Pll_qpsk_NDA(d_phi_deg, df_Rs, BlT_dB, order, EbNodB)

% il faut avoir au prealable calculer la pente du detecteur en boucle
% ouverte

load pente_NDA_QPSK
BlTT=10.^(BlT_dB);

if order==2
    zeta=sqrt(2)/2;
    A=16*zeta^2*BlTT.*(1+4*zeta^2-4*BlTT)/(1+4*zeta^2)./(1+4*zeta^2-8*zeta^2*BlTT);
    B=64*zeta^2*BlTT.^2/(1+4*zeta^2)./(1+4*zeta^2-8*zeta^2*BlTT);
elseif order==1
    B=0*BlTT;
    A=4*BlTT;
else
    display ('order1 assumed');
    B=0*BlTT;
    A=4*BlTT;
end

A=A/pente;
B=B/pente;

EbNo=10.^(EbNodB/10);
BER_th=0.5*erfc(sqrt(EbNo));

N_symb=20000;
M=4;   %QPSK
NCO_mem=0;      % initialisation NCO
filtre_mem=0;   % initialisation de la memoire du filtre
phi_est(1)=0;  % phase estimee : valeur initiale

for ii=1:N_symb

    % generation du nouveau symbole QPSK (aleatoire)

    %bits=2*randint(1,2)-1;
    bits=2*(randi(2,1,2)-1)-1;
    IE=bits(1);
    QE=bits(2);
    % symbole complexe
    symb_emis=IE+1i*QE;
    % energie du symbole
    Es=sum(abs(symb_emis).^2); % ==2 (valeur normalisee)
    
    % generation du bruit
    sigma=sqrt(Es/EbNo/4);
    noise=randn(2,1)*sigma;
    s=randn('state');
    recu=symb_emis+noise(1,:)+1i*noise(2,:); % symbole complexe bruite

    % dephasage
    d_phi=d_phi_deg*pi/180+2*pi*df_Rs*ii;   % dephasage pour le symbole no ii en radians
    recu=recu*exp(1i*d_phi);                 % symbole complexe bruite et dephase = symbole recu
    %  PLL
    out_det(ii)= -imag( (recu * exp(-1i*phi_est(ii) )).^4 );

    % filtre de boucle

    w(ii)=filtre_mem+out_det(ii); % memoire filtre + sortie detecteur
    filtre_mem=w(ii);
    out_filtre=A*out_det(ii)+B*w(ii);   % sortie du filtre a l'instant ii :  F(z)=A+B/(1-z^-1)

    %NCO
    phi_est(ii+1)=(out_filtre+NCO_mem); % N(z)=1/(z-1)
    NCO_mem=phi_est(ii+1);
end



phi_est_deg = phi_est*180/pi;
B_w = B*w/(2*pi);
% 
% figure
% plot(phi_est*180/pi)
% grid on
% xlabel('time')
% ylabel('phi-est [degre]')
% 
% figure
% plot(out_det)
% grid on
% xlabel('time')
% ylabel('detector output');
% 
% figure
% plot(B*w/(2*pi))
% grid on
% xlabel('time')
% ylabel('frequency error');