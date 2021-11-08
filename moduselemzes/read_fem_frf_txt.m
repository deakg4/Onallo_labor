function readedmatrix = read_fem_frf_txt(filename)
% filename = 'FreqResp1_DP_0.txt'; % F�jln�v megad�sa

T = readtable(filename, 'Format', 'auto'); % txt beolvas�sa t�bl�zat form�tumba
C = table2cell(T);       % t�bl�zat konvert�l�sa cella m�trixsz�
N = size(C,1);           % mintasz�m

f = C(:,1);     % frekvenciamint�k cell�kban
Amp = C(:,2);   % amplit�d�mint�k cell�kban
Phase = C(:,3); % f�zismint�k cell�kban

% a ',' karakter cser�je '.' karakterre, �s string-r�l double-re
% konvert�l�s mint�nk�nt
for i = 1:N
    
    dot_i_f = strfind(f{i},',');    % ',' karakter indexe
    f{i}(dot_i_f) = '.';            % ',' cser�je '.'-ra
    % f{i} = str2double(f{i}); % itt most nem sz�ks�ges a double-re konvert�l�s
    
    
    dot_i_Amp = strfind(Amp{i},',');
    Amp{i}(dot_i_Amp) = '.';
    Amp{i} = str2double(Amp{i}); % itt sz�ks�ges a double-re konvert�l�s

    
    dot_i_Ph = strfind(Phase{i},',');
    Phase{i}(dot_i_Ph) = '.';
    % Phase{i} = str2double(Phase{i}); % itt sem sz�ks�ges a double-re konvert�l�s
    
end

% Cell�kb�l vektorokba konvert�l�s
f = cell2mat(f);
Amp = cell2mat(Amp);
Phase = cell2mat(Phase);
readedmatrix = [f, Amp, Phase];
end

