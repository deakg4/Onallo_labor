clc
clear
close all

[x1,x2] = meshgrid(-1000:100:1000,-1000:100:1000); 
% a meshgrid functionnal lehet létrehozni olyan változókat amiket utáne 
% egy surf plothoz fel lehet használni

f =@(x1,x2) 1/7*sqrt(2500+x1.^2) + 1/4*sqrt(400+(x2-x1).^2) + 1/2*sqrt(900+(100-x2).^2);

f_2 =@(x) 1/7*sqrt(2500+x(1).^2) + 1/4*sqrt(400+(x(2)-x(1)).^2) + 1/2*sqrt(900+(100-x(2)).^2);
figure()
ezcontour(f, [-1000 1000 -1000 1000])
figure()
fy = f(x1,x2);

% surf(x1,x2,f(x1,x2))
surf(x1,x2,f(x1, x2))

fminsearch(f_2, [0 0])


