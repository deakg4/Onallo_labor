%% Clear & read
clc
clear all
close all
%addpath('M:\Documents\-----WORK-----\Matlab\LG_Gen3_housing')

% File txt-bõl való beolvasása esetén szükség van delimiterre

delimiter = 'space';
headerlinesIn = 1;
filename1_txt = 'FreqRespForce_DP_0.txt';
filename1 = 'LG_Gen3_housing_acc_Fs_20kHz_21_07_13.mat';
filename2 = 'FreqRespForce_DP_1.txt';

load(filename1)
% [fpath1, fname1, fext1] = fileparts(filename1);
% [fpath2, fname2, fext2] = fileparts(filename2);


% File beolvasása exceltáblázatból
FRF_sim_table = readtable('frf1_damped.xlsx');


FRF_sim = table2array(FRF_sim_table);
% FRF_sim = str2double(FRF_sim);

% File beolvasása txtbõl
% FRF_sim_table  = readtable(filename1_txt, 'Delimiter', delimiter, 'ReadVariableNames',false);
% FRF_sim_table2 = readtable(filename2, 'Delimiter', delimiter, 'ReadVariableNames',false);

% FRF_sim_table  = FRF_sim_table(2:5001, 1:5);
% FRF_sim_table2 = FRF_sim_table2(2:5001, 1:5);

% A másik filet is excelfileból olvastam be.
% FRF_mes_table = readtable('frf1_damped_wocube.xlsx');

j = sqrt(-1);
%% Calc

FRF_mes = FRF_matrix(1:5000,1);
FRF_mes20log = 20*log10(FRF_mes);

% Ha mind a két adat excelfileból van beolvasva
% FRF_mes = table2array(FRF_sim_table2);
% FRF_mes = str2double(FRF_mes);




% Ha txtbõl olvassuk be a filet konvertálni kell doublevé
% FRF_sim2 = table2array(FRF_sim_table2);
% FRF_sim2 = str2double(FRF_sim2);

% beolvasás excelbõl
FRF_sim_complex = FRF_sim(:,4)+(j*FRF_sim(:,5));
% FRF_mes_complex = FRF_mes(:,4)+(j*FRF_mes(:,5));
FRF_mes_complex = FRF_mes;
% beolvasás txtbõl
% FRF_sim_complex = FRF_sim(:,4)+(j*FRF_sim(:,5));


%% ---------DUMMYVECTORS---------- 
x = linspace(1,100,1000);
% dummyvector1 = 10*ones(1000,1);
% dummyvector1(495:505) = zeros(11,1);
dummyvector1 = x;
dummyvector1 = dummyvector1';
% 
% dummyvector2 = 100*ones(1000,1);
dummyvector2 = x.^2;
dummyvector2 = dummyvector2';
%  ---------DUMMYVECTORS----------

%% frac

frac = Frac(dummyvector1, dummyvector2, 1);
% frac = Frac(FRF_sim_complex ,FRF_mes_complex, 1);
% fullfrac = fullFrac(abs(FRF_sim_complex) ,abs(FRF_mes_complex), 2000, 2500);
% absfrac = abs(frac);

% frac = Frac(FRF_sim ,FRF_mes, 1);
% fullfrac = fullFrac(abs(FRF_sim) ,abs(FRF_mes), 2000, 2500);
absfrac = abs(frac);


%% Plot
close all
x = linspace(1,5000,5000);
y = linspace(1,5000,5000);

figure
% subplot(3,1,1)
% plot([20*log10(abs(FRF_mes_complex)) 20*log10(abs(FRF_sim_complex))])
% subplot(3,1,2)
% plot([((atan(imag(FRF_mes)./real(FRF_mes)))/pi*180) FRF_sim(:,3)]);

subplot(3,1,1)
plot(abs(FRF_sim_complex))
subplot(3,1,2)
plot(abs(FRF_mes_complex));

subplot(3,1,3)
plot(frac)
xlabel('Frequency [Hz]')

