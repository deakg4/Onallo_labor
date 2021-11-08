function [Uti] = elmozdulasszamitas_optimum(m, k, c, force, force_pos, Ms, omegakezdo, Nomega, Kiertekeles, visszacsatolt, i)

% Az optimumkereséshez létrehozott egyszerűbb fgv
% ez volt az eredeti fgv, de módosítom, mert nincs szükségem ennyi kimenő
% adatra belőle:
% function [U, ALFA, omega, OMEGA] = elmozdulasszamitas_optimum(m, k, c, force, force_pos, Ms, omegakezdo, Nomega, Kiertekeles, visszacsatolt)
%
%m - tömeg
%k - rugómerevség
%c - csillapítás
%force - a terhelő erő nagysáha
%Ms - tömegek száma
%force_pos - hanyadik tömeget szeretnénk terhelni
%omegakezdo - mekkora frekvenciától szeretnénk kezdeni a kiértékelést
%Nomega - a frekvenciatartomány szélessége
%Kiertekeles - az utolsó módustól vett további adatok
%visszacsatolt - a tömeg rugó rendszer visszacsatolt-e [1-igen, 0-nem]
%i- hamyas indexű elmozdulást adja visszza
%
%M - tömegmátrix
%C - csillapításmátrix
%K - rugómerevségmátrix
%FI - 1-re normált sajátvektor - Fi_n -ek módusalakok
%OMEGA - sajátértékek - omega_n-ek
%OMEGA(end)*Kiertekeles
%
%Ut - elmozdulások
%ALFA - 
%omega - frekvenciatartomány
%OMEGA - Sajátértékek
j = sqrt(-1);
F = zeros(Ms,1);
F(force_pos) = force;

[M, K, C, FI, OMEGA2] = modusmatrixgenerator(m, k, c, Ms, visszacsatolt);

KSZIOMEGAx2 = diag(FI.'*C*FI);

OMEGA = diag(sqrt(OMEGA2));
%omegakezdo = 0;
%Nomega = 1000;
%Kiertekeles = 1;
omega = linspace(omegakezdo,OMEGA(end)*Kiertekeles,Nomega).';

ALFA = zeros(Nomega,Ms);
U = zeros(Ms,Nomega);

for n = 1:Ms
    ALFA(:,n) = FI(:,n).'*F./(OMEGA(n)^2+j*omega*KSZIOMEGAx2(n)-omega.^2);
    U = U + FI(:,n).*ALFA(:,n).';
end

Ut = U';
Uti = Ut(:,i);

function [M, K, C, FI, OMEGA2] = modusmatrixgenerator(m, k, c, Ms, visszacsatolt)

    %
    %m - tömeg
    %k - rugómerevség
    %c - csillapítás
    %Ms - tömegek száma
    %visszacsatolt - a tömeg rugó rendszer visszacsatolt-e [1-igen, 0-nem]
    %M - tömegmátrix
    %C - csillapításmátrix
    %K - rugómerevségmátrix
    %FI - 1-re normált sajátvektor - Fi_n -ek módusalakok
    %OMEGA2 - sajátértékek - omega_n_negyzetek

    mm = m.*ones(1,Ms);

    cm = c*ones(1,Ms);
    km = k*ones(1,Ms);

    M = eye(Ms).*m;
    K = zeros(Ms,Ms);
    C = zeros(Ms,Ms);
    I2 = [1 -1; -1 1];

    for num = 1:Ms-1
        Ksub = km(num)*I2;
        K(num:num+1,num:num+1) = K(num:num+1,num:num+1) + Ksub;
        Csub = cm(num)*I2;
        C(num:num+1,num:num+1) = C(num:num+1,num:num+1) + Csub;
    end

    if visszacsatolt == 1
        K(Ms,1) = K(Ms,1) -km(Ms);
        K(1,Ms) = K(1,Ms) -km(Ms);
        K(1,1) = K(1,1) + km(Ms);
        K(Ms,Ms) = K(Ms,Ms) + km(Ms);
        C(Ms,1) = C(Ms,1) -cm(Ms);
        C(1,Ms) = C(1,Ms) -cm(Ms);
        C(1,1) = C(1,1) + cm(Ms);
        C(Ms,Ms) = C(Ms,Ms) + cm(Ms);
    end
        %FI - 1-re normált sajátvektor - Fi_n -ek módusalakok
        %OMEGA2 - sajátértékek - omega_n_negyzetek
        [FI,OMEGA2] = eig(M\K);
    end


end