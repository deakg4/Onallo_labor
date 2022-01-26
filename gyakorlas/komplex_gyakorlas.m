% egységgyök
clc
clear
close
i = sqrt(-1);
a = linspace(0,1,100);

b = sqrt(1 - a.^2);
z = a + b*i;

for n = 1:2000
    y(n) = z(20)^(n/10);
end

figure
h = plot(imag(y),real(y));
xlim([-1,1])
ylim([-1,1])
on = 1;
x = 0;
dir = 1;
counter = 0;

    while on == 1
        if dir == 1
            x = x+0.01;
        end
        if dir == 0
           x = x-0.01; 
        end
        b = sqrt(x - a.^2);
        z = a + b*i;
        for n = 1:2000
            y(n) = z(100)^(n/10);
        end
        set(h, 'YData', imag(y),'XData', real(y));
        pause(0.1)
        if x > 1
            dir = 0;
        end
        if x < 0
            dir = 1;
            counter = counter + 1;
        end
        if counter == 1
            on = 0;
        end
    end

    

