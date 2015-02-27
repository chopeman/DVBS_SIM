%
% This routine computes the filter parameters A and B
% for a second order loop
% The results are valid when the detector gain (G)is equal to 1
% If G~=1, the values of A and B have to be divided by G
% Note that BlT is the normalized loop bandwidth wrt loop period.

close all
clear all

BlT_dB=[-4:1:-1];

damp=[0.5,sqrt(2)/2,1];

BlT=10.^(BlT_dB);

for ii=1:length(damp)
   
zeta=damp(ii);

%wnT=2*BlT./(zeta+1/(4*zeta));
%A=wnT.*(2+zeta*wnT)./(1+3*zeta*wnT+wnT.^2);
%B=wnT.^2./(1+3*wnT+zeta*wnT.^2);

A=16*zeta^2*BlT.*(1+4*zeta^2-4*BlT)/(1+4*zeta^2)./(1+4*zeta^2-8*zeta^2*BlT);
B=64*zeta^2*BlT.^2/(1+4*zeta^2)./(1+4*zeta^2-8*zeta^2*BlT);

 
   
figure(1)

if ii==1
loglog(A,BlT,'k-')
xlabel('A')
ylabel('BlT')
title('second order loop dimensionning')
grid on
hold on
elseif ii==2
    loglog(A,BlT,'b-')
else
        loglog(A,BlT,'r-')
end

    

figure(2)
if ii==1
loglog(B,BlT,'k-')
xlabel('B')
ylabel('BlT')
title('second order loop dimensionning')
grid on
hold on
elseif ii==2
    loglog(B,BlT,'b-')
else
        loglog(B,BlT,'r-')
end

end

figure(1)
legend('zeta=0.5','zeta=sqrt(2)/2','zeta=1')

figure(2)
legend('zeta=0.5','zeta=sqrt(2)/2','zeta=1')

