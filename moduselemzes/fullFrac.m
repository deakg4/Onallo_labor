function fRac = fullFrac(H1, H2, omega1,omega2)
%FULLFRAC Summary of this function goes here
%   Detailed explanation goes here

% kiszámolja két átviteli függvény correlációját.

% omega2 - kiértékelési tartomány ha 0 akkor a H1 függvény hosszában lesz
% kiértékelve.
% H1 - elsõ átviteli függvény
% H2 - második átviteli függvény

% omega1 = 1;
if omega2 == 0
    omega2 = length(H1);
end
% A változat
    szamlalo = ((abs(H1(omega1:omega2,1)'*H2(omega1:omega2,1)))^2);
    nevezo = ((H1(omega1:omega2,1)'*H1(omega1:omega2,1))*(H2(omega1:omega2,1)'*H2(omega1:omega2,1)));
    fRac = szamlalo ./ nevezo;

end

