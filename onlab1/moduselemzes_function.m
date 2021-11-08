%moduselemzes
clc
clear
% close all
clf

%Jelmagyarázat:
%m - tömeg
%c - csillapítás
%k - rugómerevség
%M - tömegmátrix
%C - csillapításmátrix
%K - rugómerevségmátrix

%fkalap(omega) - gerjesztés
%fi_n - n-edik módushoz tartozó rezgés sajátvektor
%omega_n - módus sajátfrekvenciája
%kszi_n = c(1)/2*m*omega_n - a rendszer n-edik módusához tartozó csillapítási
%tényező

visszacsatolt = 1;
j = sqrt(-1);
%m(5) = 5; %tömeg

m = 1;
k = 1;
c = 0.00;
Ms = 18;

%f - gerjesztések
% F = zeros(Ms,1);
% F(1) = 1;
f = 1;
fnum = 1;
omegakezdo = 0;
Nomega = 10000;
Kiertekeles = 3;

[M, K, C, FI, OMEGA2] = modusmatrixgenerator(m, k, c, Ms, visszacsatolt);

[U, ALFA, omega] = elmozdulasszamitas(C, FI, OMEGA2, f, fnum, omegakezdo, Nomega, Kiertekeles, Ms);

% KSZIOMEGAx2 = diag(FI.'*C*FI);
% 
% OMEGA = diag(sqrt(OMEGA2));

% Kiertekeles = 2;
% omega = linspace(0,OMEGA(end)*Kiertekeles,Nomega).';
% 
% 
% ALFA = zeros(Nomega,Ms);
% U = zeros(Ms,Nomega);
% %F9
% 
% hold on
% for n = 1:Ms
%     ALFA(:,n) = FI(:,n).'*F./(OMEGA(n)^2+j*omega*KSZIOMEGAx2(n)-omega.^2);
%     U = U + FI(:,n).*ALFA(:,n).';
%     plot(omega,20*log10(abs(ALFA(:,n))))
% %     plot(omega,20*log10(abs(U(n,:))))
% end


hold on
for n = 1:Ms
    U = U + FI(:,n).*ALFA(:,n).';
    plot(omega,20*log10(abs(ALFA(:,n))))
%     plot(omega,20*log10(abs(U(n,:))))
end
hold off

figure
hold on
plot(omega,20*log10(abs(U(17,:))))
plot(omega,20*log10(abs(U(18,:))))
plot(omega,20*log10(abs(U(1,:))))
plot(omega,20*log10(abs(U(2,:))))
plot(omega,20*log10(abs(U(3,:))))
hold off

% hold on
% nstart = 1;
% nend = nstart+49;
% 
% for n = nstart:1:nend
%     plot(omega,20*log10(abs(U(n,:))))
% end
% hold off
