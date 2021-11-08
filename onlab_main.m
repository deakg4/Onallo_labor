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
Nomega = 10000;
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
%% Első megoldás a függvények négyzetes hibáját próbálom optimalizálni - ez
% a módszer nem vezet eredményre.

% legyen
% k = x(1)
% c = x(2)
n = 1;
p = 1;
e_real_np =@(x) (real(FRF_matrix(:,n))-real(elmozdulasszamitas_optimum(m, x(1), x(2), force, force_pos, Ms, omegakezdo, Nomega, Kiertekeles, visszacsatolt,p))).^2 ;
e_test = e_real_np([1 1]);
E_real_np =@(x) sum((real(FRF_matrix(:,n))-real(elmozdulasszamitas_optimum(m, x(1), x(2), force, force_pos, Ms, omegakezdo, Nomega, Kiertekeles, visszacsatolt,p))).^2) ;
E_test = E_real_np(1.0e+09*[1.2530 1.4251])

fminsearch(E_real_np, [3 3])
% E_real_np_to_plot =@(x1,x2) sum((real(FRF_matrix(:,n))-real(elmozdulasszamitas_optimum(m, x1, x2, force, force_pos, Ms, omegakezdo, Nomega, Kiertekeles, visszacsatolt,p))).^2) ;

% Ellenőrzés, hogy a function jó eredményeket ad e vissza
% e_test_manual = (real(FRF_matrix(:,n))-real(elmozdulasszamitas_optimum(m, 2, 3, force, force_pos, Ms, omegakezdo, Nomega, Kiertekeles, visszacsatolt,p)).^2) ;
% E_test_manual = sum((real(FRF_matrix(:,n))-real(elmozdulasszamitas_optimum(m, 2, 3, force, force_pos, Ms, omegakezdo, Nomega, Kiertekeles, visszacsatolt,p))).^2) ;
% test ha ez csak 0 ból álló vektor akkor a két megoldás megegyezik:
% e_diff = e_test_manual-e_real_np([2 3]);
% erre meg simán 0-t kell kapni:
% E_diff = E_test_manual-E_real_np([2 3]);
% A function jó eredményt ad vissza




% ki kell vonnunk egymásból a számolt elmozdulásokat és a mérési
% eredményekből kapott elmozdulásokat frekvencia tartományban.

% e_real_i = (real(FRF_matrix(:,1))-real(U(:,1))).^2;
e_real_i = (real(FRF_matrix(:,1))-real(U(:,1))).^2;
% ez egy vektort fog visszaadni ami frekvenciatartományban tartalmazza a
% hibákat.

% e_imag_i = (imag(FRF_matrix(:,1))-imag(U(:,1))).^2;
e_imag_i = (imag(FRF_matrix(:,1))-imag(U(:,1))).^2;
% ez is egy vektort ad vissza ami tartalmazza frekvenciatartományban a
% hibákat.

% többféle módszert kell felírnom
% a lényeg, hogy a két függvény különbségét kell képeznem, ezek adnak egy
% hibát, ha ezt a hibát négyzetre emelem akkor csak pozitív értékeket
% kapok. Ezeket az értékeket összeadva kaphatok egy olyan skalár
% mennyiséget ami az összes hibát leírja 
% E = sum(e_i)

% A két függvényben lehet többféle bemenet is, pl mint a két függvénynek
% veszem a valós részét képzek egy hibát, veszem külön a képzetes részt és
% annak is veszem a hibáját, ezt a két hibaértéket összeadom és ez ad egy
% teljes hibát a két függvényre.
% E_imag = sum(e_imag_i)
% E_real = sum(e_real_i)

%   pl:
% E_imag = sum(e_imag_i);
% E_real = sum(e_real_i);
% E = E_imag + E_real;

% ITT TARTOK
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% Ahhoz hogy ki tudjam rajzoltatni a bemenő paraméterek és a hiba nagysága
% között át kell alakítanom az elmozdulás számításhoz használt
% functionomet, mert jelenleg nem képes csak egyszerű változóként fogadni
% az értékeket, meg kell oldanom, hogy vagy képes legyen a vectoros
% befogadásra, vagy kell írnom egy loopot ami végigmegy az iterációkon.

% Először lehet egy looppal kezdem az egyszerűbbnek tűnik.
% Egy olyan loopra van szükségem, amiben mind a két változó végig fut egy
% egy tartományon, ehhez két iterációs változóra lesz szükségem.

%% surface generator
% legyenek az iterációs változók
% s és d
Nk = 100; % rugómerevség
Nc = 100; % csillapítási tényező
k_end = 10;
c_end = 10;


k_start = 0.01;
c_start = 0;

alteration_k = (k_end-k_start)/Nk;
alteration_c = (c_end-c_start)/Nc;

E_real_plot = zeros(Nk,Nc);
ki = 0;
ci = 0;
f = waitbar(0,'Please wait...');
for x = k_start:alteration_k:k_end
    ki=ki+1;
    ci = 0;
    for y = c_start:alteration_c:c_end
        ci=ci+1;
        E_real_plot(ki,ci) = E_real_np([x y]);        
    end
    status = (ki*ci)/(Nk*Nc);
    waitbar(status ,f ,'Please wait...');
    pause(0.01)
end
close(f)
% [X, Y] = meshgrid([0:0.1:10 0:0.1:10]);
% asd = E_real_np_to_plot(1, 1);
% surf(X, Y, E_real_np_to_plot(X, Y))
%% surface plot
X = linspace(k_start,k_end,Nk+1);
Y = linspace(c_start,c_end,Nc+1);
figure(1)
surf(X, Y, E_real_plot)



%% plots

calc_fig = figure(2);
%calc_plot = plot(abs(20*log10(abs(U(:,1)))));
calc_plot = semilogy((abs(U).^20));

figure(3)
measure_plot = semilogy((abs(FRF_matrix(:,1)).^20));

figure(4)
e_plot = plot(e_test);

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

function basicwaitbar
f = waitbar(0,'Please wait...');
pause(.5)

waitbar(.33,f,'Loading your data');
pause(1)

waitbar(.67,f,'Processing your data');
pause(1)

waitbar(1,f,'Finishing');
pause(1)

close(f)
end
%%