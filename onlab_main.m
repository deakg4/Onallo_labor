clc
clear
filename = 'FreqRespForce_DP_0.txt';
fem_frf = read_fem_frf_txt(filename);
load('LG_Gen3_housing_acc_Fs_20kHz_21_07_13.mat');

%%

%Jelmagyarázat:
%m - tömeg
%c - csillapítás
%k - rugómerevség
%M - tömegmátrix
%C - csillapításmátrix
%K - rugómerevségmátrix
%
%fkalap(omega) - gerjesztés
%fi_n - n-edik módushoz tartozó rezgés sajátvektor
%omega_n - módus sajátfrekvenciája
%kszi_n = c(1)/2*m*omega_n - a rendszer n-edik módusához tartozó csillapítási
%tényező

visual = 0; % fusson e a második programrész
visszacsatolt = 1;
j = sqrt(-1);
%m(5) = 5; %tömeg

m = 1;
k = 3;
c = 0.001;
Ms = 18;

%f - gerjesztések
% F = zeros(Ms,1);
% F(1) = 1;
f = 1;
fnum = 1;
omegakezdo = 0.1;
Nomega = 10000;
Kiertekeles = 1.5;

[M, K, C, FI, OMEGA2] = modusmatrixgenerator(m, k, c, Ms, visszacsatolt);

[U, ALFA, omega] = elmozdulasszamitas(C, FI, OMEGA2, f, fnum, omegakezdo, Nomega, Kiertekeles, Ms);
U = U';
% x(:,:,1) = FRF_matrix;
% x(:,:,2) = U;

e =@(x) (FRF_matrix(x(1),1) - U(x(2),1)).^2;
fminsearch(e, [1 1])

calc_fig = figure(1);
%calc_plot = plot(abs(20*log10(abs(U(:,1)))));
calc_plot = semilogy((abs(U(:,1))));
figure(2)
plot(abs(FRF_matrix(:,1)))
figure(3)

%%
if visual == 1
    for i = 0:0.1:10
        k = i;
        for j = 0:0.01:1
            c = j;
            [M, K, C, FI, OMEGA2] = modusmatrixgenerator(m, k, c, Ms, visszacsatolt);
            [U, ALFA, omega] = elmozdulasszamitas(C, FI, OMEGA2, f, fnum, omegakezdo, Nomega, Kiertekeles, Ms);
            U = U';
            set(calc_plot, 'YData', 20*log10(abs(U(:,1))))
            title(sprintf('k = %d, c = %d',i,j))
            pause(0.05)
        end
    end
end

%%