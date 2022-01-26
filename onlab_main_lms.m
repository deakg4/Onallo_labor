clc
clear
% close all
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

% A motor teljes tömege 282g = 0.282kg

visszacsatolt = 1; % mivel ez egy motorház egy kör alakú visszacsatolt rendszerrel lehet modellezni.
j = sqrt(-1);
mhouse = 0.282; % az LG motorház össztömege
k = 200000; % rugómerevség
c = 0; % csillapítási tényező
Ms = 18; % tömegpontok száma
m = mhouse/Ms; % a tömegpontok tömege
modal_numb = 1; % a kiértékelt tömegpont sorszáma

force = 1; % az gerjesztő erő
force_pos = 10; % a gerjesztő erő pozíciója
omegakezdo = 0.1; % a kiértékelés kezdő frekvenciája, mivel az elején 
% nagy értékről indul nem férne ki egy diagramba érdemes megvágni az elején
Nomega = 2000; % a kiértékelési frekvencia tartomány
Kiertekeles = 1.5; % a kiértékelés hossza normál esetben ezzel a szorzóval 
% szorzom a legnagyobb módust így kényelmesen belefér a plotba de most nincs rá szükségem.

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
% function [Ut] = gyorsulasszamitas_optimum(m, k, c, force, force_pos, Ms, omegakezdo, Nomega, Kiertekeles, visszacsatolt, i)
% - Ez a function számolja az amplitúdó elmozdulásokat frekvenciatartományban
% egy tömegrugó rendszerben.
% - Ennek a függvénynek a bemenő paramétereit kell optimalizálni, hogy
% közelítsük vele a mért eredményeket.
U = gyorsulasszamitas_optimum(m, k, c, force, force_pos, Ms, omegakezdo, Nomega, Kiertekeles, visszacsatolt, modal_numb);
%U_acc = (j*omega)^2;  % gyorsulás
% igazóból csak a k és c paramétert kell optimalizálnunk, ehhez a fügvény
% többi paraméterét meg kell adnunk és úgy meghívni.

%  =====================================================================
%  =====================================================================
%% Első megoldás a függvények négyzetes hibáját próbálom optimalizálni - ez
% a módszer nem vezet eredményre.

f_start = 1; % a kezdő frekvenzia [1/s]
f_end = Nomega; % a záró frekvencia [1/s]

% legyen
% k = x(1)
% c = x(2)
n = modal_numb; % a mért FRF mérési pontjának száma
p = modal_numb; % a számított FRF mérési pontjának száma

multiplier = mean(real(FRF_matrix(f_start:f_end,n))) / mean(real(U));

e_real_np =@(x) (real(FRF_matrix(f_start:f_end,n))-real(gyorsulasszamitas_optimum(m, x(1), x(2), force, force_pos, Ms, omegakezdo, Nomega, Kiertekeles, visszacsatolt,p))).^2 ;
e_test = e_real_np([1 1]);
E_real_np =@(x) sum((real(FRF_matrix(f_start:f_end,n))-real(gyorsulasszamitas_optimum(m, x(1), x(2), force, force_pos, Ms, omegakezdo, Nomega, Kiertekeles, visszacsatolt,p))).^2) ;

% E_real_np_to_plot =@(x1,x2) sum((real(FRF_matrix(:,n))-real(gyorsulasszamitas_optimum(m, x1, x2, force, force_pos, Ms, omegakezdo, Nomega, Kiertekeles, visszacsatolt,p))).^2) ;

% Ellenőrzés, hogy a function jó eredményeket ad e vissza
% e_test_manual = (real(FRF_matrix(:,n))-real(gyorsulasszamitas_optimum(m, 2, 3, force, force_pos, Ms, omegakezdo, Nomega, Kiertekeles, visszacsatolt,p)).^2) ;
% E_test_manual = sum((real(FRF_matrix(:,n))-real(gyorsulasszamitas_optimum(m, 2, 3, force, force_pos, Ms, omegakezdo, Nomega, Kiertekeles, visszacsatolt,p))).^2) ;
% test ha ez csak 0 ból álló vektor akkor a két megoldás megegyezik:
% e_diff = e_test_manual-e_real_np([2 3]);
% erre meg simán 0-t kell kapni:
% E_diff = E_test_manual-E_real_np([2 3]);
% A function jó eredményt ad vissza

% ki kell vonnunk egymásból a számolt elmozdulásokat és a mérési
% eredményekből kapott elmozdulásokat frekvencia tartományban.

% e_real_i = (real(FRF_matrix(:,1))-real(U(:,1))).^2;
e_real_i = (real(FRF_matrix(f_start:f_end,1))-real(U(:,1))).^2;
% ez egy vektort fog visszaadni ami frekvenciatartományban tartalmazza a
% hibákat.

% e_imag_i = (imag(FRF_matrix(:,1))-imag(U(:,1))).^2;
e_imag_i = (imag(FRF_matrix(f_start:f_end,1))-imag(U(:,1))).^2;
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
% ki és ci
Nk = 40; % rugómerevség
Nc = 40; % csillapítási tényező

k_start = 0; % k kezdőértéke
k_end = 10000000; % k utolsó értéke

c_start = 0; % c kezdőértéke
c_end = 800; % c utolsó értéke


alteration_k = (k_end-k_start)/Nk; % k osztása
alteration_c = (c_end-c_start)/Nc; % c osztása

surface_plot = 1;
if surface_plot == 1
    E_real_plot = zeros(Nk,Nc); % egy Nk*Nc méretű mátrix a kapott eredményeknek.
    ki = 0; % k futóváltozója
    ci = 0; % c futóváltozója
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
end

%%
% a kapott eredményt érdemes lenne visszahelyettesíteni az elmozdulás
% számító functionbe, mert akkor vizuálisan is lehetne látni a kiadott
% minimumpontokon felvett érték hogyan néz ki.

f_min = fminsearch(E_real_np, [k c])
E_min = E_real_np(f_min)

U_optimum = gyorsulasszamitas_optimum(m, f_min(1), f_min(2), force, force_pos, Ms, omegakezdo, Nomega, Kiertekeles, visszacsatolt,modal_numb);

%% surface plot
X = linspace(k_start,k_end,Nk+1);
Y = linspace(c_start,c_end,Nc+1);
figure(1)
surf(X, Y, 20*log10(E_real_plot))
title('error surface')
xlabel('k')
ylabel('c')
zlabel('error')

%% plots

figure(3);
calc_plot = plot(20*log10(abs(U)));
title('Calculated frf')

figure(4)
measure_plot = plot(20*log10(abs(FRF_matrix(f_start:f_end,1))));
title('Measured frf')

figure(5)
e_real_i_plot = plot(20*log10(e_real_i));
title('e real')

figure(6)
optimum_plot_1 = plot(20*log10(abs(U_optimum)));
title('Optimum')

% figure(7)
% optimum_plot_2 = semilogy((abs(U_optimum_2).^20));
% title('frac product')

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

% További ötletek:
% Meg kell keresni az frf lokális maximumait ezek a lokális maximumok
% ugyebár a sajátmódusok helyein lesznek mert ott száll el a gyorsulásunk.
% Erre kell keresni valami fasza kis megoldást, hogyan lehet lokális
% amximumokat keresni
% FELADAT:
% Lokális maximum keresésének irodalmazása.
