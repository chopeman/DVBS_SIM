function [pente] = PLL_QPSK_BO_NDA(d_phi_deg, EbNodB)

EbNo=10.^(EbNodB/10);

N_symb=1000;

for jj=1:length(d_phi_deg)  % boucle sur erreur de phase
    % affichage phases
    if mod(jj,5)==0
        d_phi_deg(jj);
    end
    
    for ii=1:N_symb   % boucle sur symboles
        %bits=2*randint(1,2)-1;
        bits=2*(randi(2,1,2)-1)-1;
        IE=bits(1);
        QE=bits(2);
        symb_emis=IE+1i*QE;
        Es=sum(abs(symb_emis).^2);
        %bruit
        sigma=sqrt(Es/EbNo/4);
        noise=randn(2,1)*sigma;
        %rajouter l'erreur de phase + le bruit
        recu = symb_emis.*exp(1i*deg2rad(d_phi_deg(jj))) + noise(1) + 1i*noise(2);
        %  detecteur :
        out_det(ii) = -imag(recu.^4);
    end
    
    S_curve(jj)=mean(out_det);

end

% on calcule la pente de la caracteristique (S-Curve) autour de 0 (entre -3 et 3 degres) :
pente=S_curve((length(S_curve)+1)/2+3)-S_curve((length(S_curve)+1)/2-3);
pente=pente/(6*(d_phi_deg(2)-d_phi_deg(1))*pi/180);


figure
plot(d_phi_deg,S_curve,'b-')
grid on
hold on
plot([-3:1:3],[-3:1:3]*pente*pi/180,'r-');
xlabel('Phase error')
ylabel('Detector output')
title('Characteristic (S-curve) of detector')

save pente_NDA_QPSK pente




