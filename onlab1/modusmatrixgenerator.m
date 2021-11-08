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

mm = m*ones(1,Ms);

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