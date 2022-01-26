function fRac = fullFrac(H1, H2, omega1,omega2)
%FULLFRAC calculate correlation bethween two frequency response functions.
%   Detailed explanation goes here
%
% omega1 - ki�rt�kel�si tartom�ny kezdete
% omega2 - ki�rt�kel�si tartom�ny ha 0 akkor a H1 f�ggv�ny hossz�ban lesz
% ki�rt�kelve.
% H1 - els� �tviteli f�ggv�ny
% H2 - m�sodik �tviteli f�ggv�ny

 if ~exist('omega1','var')
     % third parameter does not exist, so default it to something
      omega1 = 1;
 end

 if ~exist('omega2','var')
     % third parameter does not exist, so default it to something
      omega2 = length(H1);
 end

% A v�ltozat
    szamlalo = ((abs(H1(omega1:omega2,1)'*H2(omega1:omega2,1)))^2);
    nevezo = ((H1(omega1:omega2,1)'*H1(omega1:omega2,1))*(H2(omega1:omega2,1)'*H2(omega1:omega2,1)));
    fRac = szamlalo ./ nevezo;

end

