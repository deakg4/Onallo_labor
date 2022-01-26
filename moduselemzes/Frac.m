function fRac = Frac(H1, H2, ablak)
% This function calculates frac between H1 and H2
%
% szamlalo = (abs(H1*H2'))^2;
% nevezo = (H1*H1')*(H2*H2');

omega1 = 1;
omega2 = length(H1);
fRac = zeros(omega2-ablak,1);
for omega = omega1:(omega2-ablak)

    
% B változat
    
    szamlalo = ((abs(H1(omega:omega+ablak,1)'*H2(omega:omega+ablak,1)))^2);
    nevezo = ((H1(omega:omega+ablak,1)'*H1(omega:omega+ablak,1))*(H2(omega:omega+ablak,1)'*H2(omega:omega+ablak,1)));
    fRac(omega) = szamlalo / nevezo;
%     fRac(omega) = ((abs(H1(omega1:omega,1)* H2(omega1:omega,1)'))^2) / ((H1(omega1:omega,1)*H1(omega1:omega,1)')*(H2(omega1:omega,1)*H2(omega1:omega,1)'));
end
% fRac = ((abs(H1*H2'))^2)./((H1*H1')*(H2*H2'));

end

