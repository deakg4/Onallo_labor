function fRac = fullFrac(H1, H2, omega1,omega2)
%FULLFRAC Summary of this function goes here
%   Detailed explanation goes here

% kisz�molja k�t �tviteli f�ggv�ny correl�ci�j�t.

% omega2 - ki�rt�kel�si tartom�ny ha 0 akkor a H1 f�ggv�ny hossz�ban lesz
% ki�rt�kelve.
% H1 - els� �tviteli f�ggv�ny
% H2 - m�sodik �tviteli f�ggv�ny

% omega1 = 1;
if omega2 == 0
    omega2 = length(H1);
end
% A v�ltozat
    szamlalo = ((abs(H1(omega1:omega2,1)'*H2(omega1:omega2,1)))^2);
    nevezo = ((H1(omega1:omega2,1)'*H1(omega1:omega2,1))*(H2(omega1:omega2,1)'*H2(omega1:omega2,1)));
    fRac = szamlalo ./ nevezo;

end

