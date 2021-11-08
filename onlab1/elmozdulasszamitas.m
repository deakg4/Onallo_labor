function [U, ALFA, omega, OMEGA] = elmozdulasszamitas(C, FI, OMEGA2, f, fnum, omegakezdo, Nomega, Kiertekeles, Ms)

j = sqrt(-1);
F = zeros(Ms,1);
F(fnum) = f;

KSZIOMEGAx2 = diag(FI.'*C*FI);

OMEGA = diag(sqrt(OMEGA2));
%omegakezdo = 0;
%Nomega = 1000;
%Kiertekeles = 1;
omega = linspace(omegakezdo,OMEGA(end)*Kiertekeles,Nomega).';

ALFA = zeros(Nomega,Ms);
U = zeros(Ms,Nomega);

for n = 1:Ms
    ALFA(:,n) = FI(:,n).'*F./(OMEGA(n)^2+j*omega*KSZIOMEGAx2(n)-omega.^2);
    U = U + FI(:,n).*ALFA(:,n).';
end

end