function readedmatrix = read_fem_frf_txt(filename)
% filename = 'FreqResp1_DP_0.txt'; % Fájlnév megadása

T = readtable(filename, 'Format', 'auto'); % txt beolvasása táblázat formátumba
C = table2cell(T);       % táblázat konvertálása cella mátrixszá
N = size(C,1);           % mintaszám

f = C(:,1);     % frekvenciaminták cellákban
Amp = C(:,2);   % amplitúdóminták cellákban
Phase = C(:,3); % fázisminták cellákban

% a ',' karakter cseréje '.' karakterre, és string-rõl double-re
% konvertálás mintánként
for i = 1:N
    
    dot_i_f = strfind(f{i},',');    % ',' karakter indexe
    f{i}(dot_i_f) = '.';            % ',' cseréje '.'-ra
    % f{i} = str2double(f{i}); % itt most nem szükséges a double-re konvertálás
    
    
    dot_i_Amp = strfind(Amp{i},',');
    Amp{i}(dot_i_Amp) = '.';
    Amp{i} = str2double(Amp{i}); % itt szükséges a double-re konvertálás

    
    dot_i_Ph = strfind(Phase{i},',');
    Phase{i}(dot_i_Ph) = '.';
    % Phase{i} = str2double(Phase{i}); % itt sem szükséges a double-re konvertálás
    
end

% Cellákból vektorokba konvertálás
f = cell2mat(f);
Amp = cell2mat(Amp);
Phase = cell2mat(Phase);
readedmatrix = [f, Amp, Phase];
end

