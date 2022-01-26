function fRac = fullFrac(H1, H2, omega1,omega2)
%FULLFRAC calculate correlation bethween two frequency response functions.
%   Detailed explanation goes here
%
% omega1 - kiértékelési tartomány kezdete
% omega2 - kiértékelési tartomány ha 0 akkor a H1 függvény hosszában lesz
% kiértékelve.
% H1 - elsõ átviteli függvény
% H2 - második átviteli függvény

 if ~exist('omega1','var')
     % third parameter does not exist, so default it to something
      omega1 = 1;
 end

 if ~exist('omega2','var')
     % third parameter does not exist, so default it to something
      omega2 = length(H1);
 end

% A változat
    szamlalo = ((abs(H1(omega1:omega2,1)'*H2(omega1:omega2,1)))^2);
    nevezo = ((H1(omega1:omega2,1)'*H1(omega1:omega2,1))*(H2(omega1:omega2,1)'*H2(omega1:omega2,1)));
    fRac = szamlalo ./ nevezo;

end

