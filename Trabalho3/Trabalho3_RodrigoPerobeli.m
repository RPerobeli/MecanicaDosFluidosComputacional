clear all
close all
clc

%% Condições iniciais e de contorno;

xi=0;
xf=3;
numElem = 100;
if(rem(numElem,2)~= 0)
    numElem = numElem +1;
end
numNos = numElem+1;
dx = (xf-xi)/numElem;
x = [xi:dx:xf];
A = 1+2.2*(x-1.5).^2;
Courant = 1.2;
gamma = 1.4;

%cond iniciais
rho = 1-0.3146.*x;
T = 1-0.2314.*x;
for(j=1:numNos)
    V(j) = (0.1+1.09.*x(j))*T(j).^(0.5);
end

rho_barra = zeros(1,numNos);
V_barra = zeros(1,numNos);
T_barra = zeros(1,numNos);

critParada=1;

dRho_p = zeros(1,numNos);
dV_p = zeros(1,numNos);
dT_p = zeros(1,numNos);

dRho_c = zeros(1,numNos);
dV_c = zeros(1,numNos);
dT_c = zeros(1,numNos);

cont =1;
while(critParada > 1e-6)
    dt_vetor = Courant*dx*(T.^0.5+V).^-1; % calcular a cada iteração
    dt = min(dt_vetor);
    %% Passo Preditor
   
    %condiçoes iniciais de contorno
    rho(1) = 1;
    V(1) = 2*V(2)-V(3);
    T(1) = 1;

    rho(numNos) = 2*rho(numNos-1)- rho(numNos-2);
    V(numNos) = 2*V(numNos-1) - V(numNos-2);
    T(numNos) = 2*T(numNos-1) - T(numNos-2);
    for (i=2:numNos-1)
        %preditor da densidade
        dRho_p(i) = (-rho(i)/dx)*(V(i+1)-V(i)) - (((rho(i)*V(i))/dx)*(log(A(i+1))-log(A(i)))) - (V(i)/dx)*(rho(i+1)-rho(i));
        rho_barra(i) = rho(i) + dRho_p(i)*dt;
        
        %preditor da velocidade
        dV_p(i) = -(V(i)/dx)*(V(i+1)-V(i)) - (1/gamma)*((T(i+1)-T(i))/dx + (T(i)/(rho(i)*dx))*(rho(i+1)-rho(i)));
        V_barra(i) = V(i)+ dt*dV_p(i);
        %preditor da temperatura
        dT_p(i) = -(V(i)/dx)*(T(i+1)-T(i)) - ((gamma-1)*T(i))/dx*((V(i+1)-V(i)) + V(i)*(log(A(i+1))-log(A(i))));
        T_barra(i) = T(i) + dT_p(i)*dt;
    end
    
    %% Passo Corretor
    
    rho_barra(1) = 1;
    V_barra(1) = 2*V_barra(2)-V_barra(3);
    T_barra(1) = 1;

    rho_barra(numNos) = 2*rho_barra(numNos-1)- rho_barra(numNos-2);
    V_barra(numNos) = 2*V_barra(numNos-1) - V_barra(numNos-2);
    T_barra(numNos) = 2*T_barra(numNos-1) - T_barra(numNos-2);
    
    for (i=2:numNos-1)
        %preditor da densidade
        dRho_c(i) = (-rho_barra(i)/dx)*(V_barra(i)-V_barra(i-1)) - (((rho_barra(i)*V_barra(i))/dx)*(log(A(i))-log(A(i-1)))) - (V_barra(i)/dx)*(rho_barra(i)-rho_barra(i-1));

        %preditor da velocidade
        dV_c(i) = -(V_barra(i)/dx)*(V_barra(i)-V_barra(i-1)) - (1/gamma)*((T_barra(i)-T_barra(i-1))/dx + (T_barra(i)/(rho_barra(i)*dx))*(rho_barra(i)-rho_barra(i-1)));
        
        %preditor da temperatura
        dT_c(i) = -(V_barra(i)/dx)*(T_barra(i)-T_barra(i-1)) - ((gamma-1)*T_barra(i))/dx*((V_barra(i)-V_barra(i-1)) + V_barra(i)*(log(A(i))-log(A(i-1))));
    end
    
    dRho = 0.5*(dRho_p + dRho_c);
    dV = 0.5*(dV_p + dV_c);
    dT = 0.5*(dT_p + dT_c);
    
    rho = rho + dRho*dt;
    V = V+ dV*dt;
    T = T+ dT*dt;
    
    %Armazenamento do bocal no tempo:
    rho_bocal(cont) = rho(numElem/2+1);
    V_bocal(cont) = V(numElem/2+1);
    T_bocal(cont) = T(numElem/2+1);
    cont = cont+1;
    
    erro = [max(abs(dRho)),max(abs(dV)), max(abs(dT))];
    %critParada = critParada - 0.5;
    critParada(cont) = max(erro);
end

%% Valores no Bocal

rho_bocal_exato = (1+(gamma-1)/2)^(-1/(gamma-1));
T_bocal_exato = (1+(gamma-1)/2)^(-1);

erroRho = (abs(rho_bocal_exato - rho_bocal(cont-1)))/rho_bocal_exato*100
erroT = (abs(T_bocal_exato - T_bocal(cont-1)))/T_bocal_exato*100

erroV = (abs(0.930 - V_bocal(cont-1)))/0.930*100

%% Plots

figure;
plot(x,rho, 'b');
xlabel("X");
ylabel("\rho");
title("Densidade ao longo do bocal");

figure;
plot(x,V, 'r');
xlabel("X");
ylabel("V");
title("Velocidade ao longo do bocal");

figure;
plot(x,T, 'g');
xlabel("X");
ylabel("T");
title("Temperatura ao longo do bocal");

figure;
plot(1:cont-1,rho_bocal, 'b');
xlabel("iterações");
ylabel("\rho");
title("\rho no Bocal ao longo do tempo");

figure;
plot(1:cont-1,V_bocal, 'r');
xlabel("iterações");
ylabel("V");
title("V no Bocal ao longo do tempo");

figure;
plot(1:cont-1,T_bocal, 'g');
xlabel("iterações");
ylabel("T");
title("T no Bocal ao longo do tempo");

figure;
plot(1:cont,critParada);
xlabel("iterações");
ylabel("erro");
title("Critério de parada");
