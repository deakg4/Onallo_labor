clc
clear all
x = linspace(-10, 10, 1000)
y = exp(1i*x)
z = exp(-1i*x)
yz = y + z
plot3(x, real(y), imag(y))
hold on
plot3(x, real(z), imag(z))
plot3(x, real(yz), imag(yz))
hold off

legend({'y', 'z', 'yz'})
x_abs = abs(x)
y_abs = abs(y)