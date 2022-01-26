clc
clear
close all
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

visszacsatolt = 1;
j = sqrt(-1);

m = 1;
% 1.0e+09*[1.2530 1.4251]
k = 1;
c = 1;
Ms = 18;

%f - gerjesztések
% F = zeros(Ms,1);
% F(1) = 1;
force = 1;
force_pos = 1;
omegakezdo = 0.1;
Nomega = 5000;
Kiertekeles = 1.5;

% legenerálja a megadott paraméterekkel a mátrixokat
% ez volt a régifgv de lecserélem egy másikra
% [M, K, C, FI, OMEGA2] = modusmatrixgenerator(m, k, c, Ms, visszacsatolt);

% kiszámolja a sajátértékeket és sajátvektorokat és az ezlmozdulásokat
% frekvenciatartományban.
% ez volt a régi fgv de lecserélem egy másikra
% [U, ALFA, omega] = elmozdulasszamitas(C, FI, OMEGA2, f, fnum, omegakezdo, Nomega, Kiertekeles, Ms);

% a régi fgv-ekkel szükség van egy transponálással mert az máshogy adja
% vissza az eredményeket:
%U = U';

% - új fgv, csak az elmozdulásokat adja vissza és egybevontam a két fgv:
% function [Ut] = elmozdulasszamitas_optimum(m, k, c, force, force_pos, Ms, omegakezdo, Nomega, Kiertekeles, visszacsatolt, i)
% - Ez a function számolja az amplitúdó elmozdulásokat frekvenciatartományban
% egy tömegrugó rendszerben.
% - Ennek a függvénynek a bemenő paramétereit kell optimalizálni, hogy
% közelítsük vele a mért eredményeket.
U = elmozdulasszamitas_optimum(m, k, c, force, force_pos, Ms, omegakezdo, Nomega, Kiertekeles, visszacsatolt,3);
% igazóból csak a k és c paramétert kell optimalizálnunk, ehhez a fügvény
% többi paraméterét meg kell adnunk és úgy meghívni.

%  =====================================================================
%  =====================================================================

%% Második megoldás FRAC
% a FRAC - Frequency Assurance Criterion function a function whichi is
% compare two different frequency respons or any other type functions.

% ezzel a módszerrel a két függvényt teljes frekvencia tartományban
% vizsgálom és a két vektor által bezárt szög cosinusát elemzem, ha az
% eredmény 1-azt jelenti cos(0) = 1 tehát a két függvény párhuzamos
% egymással. Ebben az esetben 100% a korreláció
% Ha az eredmény 0 a két függvény egymással 90 fokos szöget zár be a
% korreláció 0.
% az x vector fogja tartalmazni a k és c bemenő paramétereket.
f_start = 1;
f_end = 5000;
fullfrac = fullFrac(FRF_matrix(f_start:f_end,1), elmozdulasszamitas_optimum(m, k, c, force, force_pos, Ms, omegakezdo, Nomega, Kiertekeles, visszacsatolt,3));
% ez így egy 0-hoz közeli értéket adott vissza.
% ezt a fullfrac függvényt kell ábrázolnom k és c paraméterek mentén és az
% így kapott felület írja le a hibafelületet.
frac = Frac(FRF_matrix(f_start:f_end,1), elmozdulasszamitas_optimum(m, k, c, force, force_pos, Ms, omegakezdo, Nomega, Kiertekeles, visszacsatolt,3),3);
mean_frac = mean(frac);
prod_frac = prod(frac);

E_mean_frac_np =@(x) mean(Frac(FRF_matrix(f_start:f_end,1), elmozdulasszamitas_optimum(m, x(1), x(2), force, force_pos, Ms, omegakezdo, Nomega, Kiertekeles, visszacsatolt,3),3));
E_prod_frac_np = @(x) prod(Frac(FRF_matrix(f_start:f_end,1), elmozdulasszamitas_optimum(m, x(1), x(2), force, force_pos, Ms, omegakezdo, Nomega, Kiertekeles, visszacsatolt,3),3));


%% surface generator
% legyenek az iterációs változók
% s és d
Nk = 20; % rugómerevség
Nc = 20; % csillapítási tényező
k_end = 200;
c_end = 2;

k_start = 0.1;
c_start = 0.01;

alteration_k = (k_end-k_start)/(Nk-1);
alteration_c = (c_end-c_start)/(Nc-1);

E_mean_frac_plot = zeros(Nk,Nc);
ki = 0;
ci = 0;
f = waitbar(0,'Calculate...');
for x = k_start:alteration_k:k_end
    ki=ki+1;
    ci = 0;
    for y = c_start:alteration_c:c_end
        ci=ci+1;
        E_mean_frac_plot(ki,ci) = E_mean_frac_np([x y]);        
    end
    status = (ki*ci)/(Nk*Nc);
    waitbar(status ,f ,'Calculate...');
    pause(0.01)
end
close(f)


E_prod_frac_plot = zeros(Nk,Nc);
ki = 0;
ci = 0;
f = waitbar(0,'Calculate...');
for x = k_start:alteration_k:k_end
    ki=ki+1;
    ci = 0;
    for y = c_start:alteration_c:c_end
        ci=ci+1;
        E_prod_frac_plot(ki,ci) = E_prod_frac_np([x y]);        
    end
    status = (ki*ci)/(Nk*Nc);
    waitbar(status ,f ,'Calculate...');
    pause(0.01)
end
close(f)

%% surface plot

X = linspace(k_start,k_end-k_start,Nk);
Y = linspace(c_start,c_end-c_start,Nc);
% figure(1)
% surf_frac_prod = surf(X, Y, E_prod_frac_plot);
% ax1 = gca
% title(ax1, 'Frac product surface')

figure(2)
surf_frac_mean = surf(X, Y, E_mean_frac_plot);
title('Frac mean surface')

%%
% a kapott eredményt érdemes lenne visszahelyettesíteni az elmozdulás
% számító functionbe, mert akkor vizuálisan is lehetne látni a kiadott
% minimumpontokon felvett érték hogyan néz ki.

f_mean_min = fminsearch(E_mean_frac_np, [200 0])
e_mean_min = E_mean_frac_np(f_mean_min)
f_prod_min = fminsearch(E_prod_frac_np, [200 0])
e_prod_min = E_prod_frac_np(f_prod_min)

U_optimum_1 = elmozdulasszamitas_optimum(m, f_mean_min(1), f_mean_min(2), force, force_pos, Ms, omegakezdo, Nomega, Kiertekeles, visszacsatolt,3);
U_optimum_2 = elmozdulasszamitas_optimum(m, f_prod_min(1), f_prod_min(2), force, force_pos, Ms, omegakezdo, Nomega, Kiertekeles, visszacsatolt,3);




%% plots

figure(3);
calc_plot = semilogy((abs(U).^20));
title('Calculated frf')

figure(4)
measure_plot = plot(20*log10(abs(FRF_matrix(f_start:f_end,1))));
title('Measured frf')

figure(5)
frac_plot = plot(frac);
title('frac')

figure(6)
optimum_plot_1 = semilogy((abs(U_optimum_1).^20));
title('frac mean')

figure(7)
optimum_plot_2 = semilogy((abs(U_optimum_2).^20));
title('frac product')

%%
visual = 0;
if visual == 1
    for i = 0:0.1:10
        k = i;
        for j = 0:0.01:1
            c = j;
            [M, K, C, FI, OMEGA2] = modusmatrixgenerator(m, k, c, Ms, visszacsatolt);
            [U, ALFA, omega] = elmozdulasszamitas(C, FI, OMEGA2, force, force_pos, omegakezdo, Nomega, Kiertekeles, Ms);
            U = U';
            set(calc_plot, 'YData', 20*log10(abs(U(:,1))))
            title(sprintf('k = %d, c = %d',i,j))
            pause(0.05)
        end
    end
end
