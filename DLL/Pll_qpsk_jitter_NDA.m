close all

clear all



load pente_NDA_QPSK



BlT_dB=input('log10(BlT)=[..] (ex: -3 pour 0.1%)?');  % ex: -3

BlTT=10.^(BlT_dB);

ordre=input('loop order (1 or 2)?'); 



if ordre==2

zeta=sqrt(2)/2;

A=16*zeta^2*BlTT.*(1+4*zeta^2-4*BlTT)/(1+4*zeta^2)./(1+4*zeta^2-8*zeta^2*BlTT);

B=64*zeta^2*BlTT.^2/(1+4*zeta^2)./(1+4*zeta^2-8*zeta^2*BlTT);



elseif ordre==1

   B=0*BlTT;

   A=4*BlTT;

else

   display ('order1 assumed');

    B=0*BlTT;

    A=4*BlTT;

end

 

 A=A/pente;

 B=B/pente;

 

EbNodB=input('Eb/No dB=?');

EbNo=10.^(EbNodB/10);



N_symb=1e5;

M=4;   %QPSK





for jj=1:length(BlTT) % loop on BlT

   jj

   BlT=BlTT(jj);

NCO_mem=0;

filtre_mem=0;

phi_est(1)=0;



for ii=1:N_symb



    if mod(ii,1000)==0

        ii

    end

    

    

   %bits=2*randint(1,2)-1;

   bits=2*(randi(2,1,2)-1)-1;

   IE=bits(1);

   QE=bits(2);

   symb_emis=IE+j*QE;

 

 

 Es=sum(abs(symb_emis).^2);

 

 

 %

 %bruit

 %

   sigma=sqrt(Es/EbNo/4);

   

   noise=randn(2,1)*sigma;

 

   recu=symb_emis+noise(1,:)+j*noise(2,:);

   

   %

   % dephasage

   %

   

   d_phi=0;

   recu=recu*exp(j*d_phi);

   

   %

   %  PLL

   %

   

   out_det(ii)=-imag((recu*exp(-j*phi_est(ii)))^4);

   

   % filtre de boucle

   

   w(ii)=filtre_mem+out_det(ii);

   filtre_mem=w(ii);

   out_filtre=A(jj)*out_det(ii)+B(jj)*w(ii);

   

   %NCO

   

   phi_est(ii+1)=out_filtre+NCO_mem;

   NCO_mem=phi_est(ii+1);

 

  

end



error=phi_est(5000:length(phi_est));

jitter(jj)=var(error);



end



% plotting of results



figure(1)

loglog(BlTT,jitter,'k*-')

xlabel('BlT')

ylabel('jitter [rad2]')

title('jitter ')

grid on



figure(2)

loglog(BlTT,sqrt(jitter)*180/pi,'k*-')

grid on

xlabel('BlT')

ylabel('std deviation [deg]')

title('ecart type en degre')



figure(3)

semilogx(BlTT,sqrt(jitter)*180/pi,'k*-')

grid on

xlabel('BlT')

ylabel('std deviation [deg]')

title('ecart type en degre')




