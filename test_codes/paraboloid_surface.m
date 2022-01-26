% paraboloid

f_paraboloid  =@(x) -((x(1)-0)^2 + (x(2)-0)^2) + 4;

f_paraboloid_inv =@(x) (-1)*(f_paraboloid(x)); 

N = 100;
x_start = -10;
x_end = 10;
y_start = -10;
y_end = 10;

X = linspace(x_start, x_end, N);
Y = linspace(y_start, y_end, N);

x_alternation = (x_end-x_start)/(N-1);
y_alternation = (y_end-y_start)/(N-1);
xi = 0;
yi = 0;

surf_paraboloid = zeros(N,N);

for x = x_start:x_alternation:x_end
    xi = xi+1;
    yi = 0;
    for y = y_start:y_alternation:y_end
        yi = yi+1;
        surf_paraboloid(xi, yi) = f_paraboloid([x, y]);
        surf_paraboloid_inv(xi, yi) = f_paraboloid_inv([x, y]);
    end
end

figure(1)
surf(X, Y, surf_paraboloid)
figure(2)
surf(X, Y, surf_paraboloid_inv)
paraboloidmax = fminsearch(f_paraboloid_inv, [1 1])

