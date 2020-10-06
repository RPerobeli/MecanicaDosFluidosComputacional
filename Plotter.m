close all
clc
clear all

numElem = 100;
numNos = numElem+1;
xi = -2;
xf = 2;
dx = (xf-xi)/numElem;

x = [xi:dx:xf];

for i=1:numNos
    if(x(i)<-1 && x(i)>=-2)
        y(i) = 1;
        z(i) = -1;
    end
    if(x(i)<0 && x(i)>=-1)
        y(i) = 0;
        z(i) = 0;
    end
    if(x(i)<1 && x(i)>=0)
        y(i) = 0;
        z(i) = 0;
    end
    if(x(i)<=2 && x(i)>=1)
        y(i) = 1;
        z(i) = 1;
    end
end

figure;
plot(x,y,'b');
xlabel("x");
ylabel("f(x)");
title("Prolongamento Par");

figure
plot(x,z,'r');
xlabel("x");
ylabel("f(x)");
title("Prolongamento Ímpar");
