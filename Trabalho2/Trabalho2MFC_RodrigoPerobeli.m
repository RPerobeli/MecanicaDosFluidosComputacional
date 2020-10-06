% Inicio
clear all
close all
clc

%Dominio temporal:
tempoInicial = 0.0;             %MUDAR
%tempoFinal = 0.5;
tempoFinal = 1.0;
%Dominio espacial
Xi = -1.0;
Xf = 1.0;                       %MUDAR
Yi = -1.0;
Yf = 1.0;
%Definição do número de repetições para convergência (4,8,16,32,64)
NumMalhas = 5; 

 for k = 1:NumMalhas
    %O número de elementos deve ser definido de forma que possa ser variável
    %para depois fazer a convergência
    numElementos = 2^(k+1);
%     numElementos = 32;
    numNos = numElementos+1;
    %Dimensão do elemento
    h = (Xf-Xi)/numElementos;
    h2 = h*h;
    %passo temporal
    dt = h2;
    %numero de Reynolds
    Re = 20;                       %MUDAR
    %Re =1;
    
    %coordenadas:
    x = zeros(numNos,1);
    y = zeros(numNos,1);
    x(1) = Xi;
    y(1) = Yi;
    for i = 2:numNos
        %Define um vetor de coordenadas 
        x(i) = Xi+(i-1)*h;
        y(i) = Yi+(i-1)*h;
    end
    %Condição inicial
    t= tempoInicial;
    for i = 1:numNos
        for j = 1:numNos
            %Solução exata (abrange inicial + contorno)
            %w(i,j) = sin(pi*x(i))*sin(pi*y(j));
            w(i,j) = 2*pi*cos(pi*x(i))*cos(pi*y(j));           %MUDAR
        end
    end
    %sigma é um numero de CFL,
    Pe=1/Re; %Pe = epsilon = numero de pecle
    sigma = Pe*(dt/h2);
    
    while t < tempoFinal
        %resolve para cada tempo a malha inteira
        t = t+dt/2.;
        
        %velocidade de convecção
        u1 = (-cos(pi*x(:))*sin(pi*y(:))')*exp(-2*pi*pi*t/Re);                  %MUDAR
        u2 = (sin(pi*x(:))*cos(pi*y(:))')*exp(-2*pi*pi*t/Re);
        %u1 = zeros(numNos,numNos);
        %u2 = zeros(numNos,numNos);
        
        rho1 = u1*dt/h;
        rho2 = u2*dt/h;
        
        % condicoes de contorno
        for cont = 1:numNos
            %X:
            w_meio(1,cont) = 2*pi*cos(pi*x(1))*cos(pi*y(cont))*exp(-2*pi*pi*t/Re);
            w_meio(numNos,cont) = 2*pi*cos(pi*x(numNos))*cos(pi*y(cont))*exp(-2*pi*pi*t/Re);                
            %Y:
            w_meio(cont,1) = 2*pi*cos(pi*x(cont))*cos(pi*y(1))*exp(-2*pi*pi*t/Re);
            w_meio(cont,numNos) = 2*pi*cos(pi*x(cont))*cos(pi*y(numNos))*exp(-2*pi*pi*t/Re);
        end
        
        
        %METODO ADI
        %primeiro passo do ADI implicito em X
        for j = 2:numNos-1
            for i = 2:numNos-1
                a(1,i-1) = -1*(sigma/2+(rho1(i,j)/4));%diagonal inferior
                b(1,i-1) = (1+sigma); %diagonal principal
                c(1,i-1) = ((-sigma/2)+(rho1(i,j)/4));%diagonal superior           
                %notar que o rho1 e rho2 variam de 2 até numNos-1
                f(i-1)= (1-sigma)*w(i,j)+ (sigma/2-rho2(i,j)/4)*w(i,j+1)+ ((sigma/2)+(rho2(i,j)/4))*w(i,j-1);
            end
            %COLOCAR AQUI OS TERMOS FONTE QDO i=1 E i=4
            f(1) = f(1)+ w_meio(1,j)*((sigma/2)+(rho1(2,j)/4)); %rho1(1,j) faz parte do contorno
            f(numNos-2) = f(numNos-2) - ((-sigma/2)+(rho1(numNos-1,j)/4))*w_meio(numNos,j); %rho1(numNos,j) faz parte do contorno
            w_meio(2:numNos-1,j) = TDMAsolver(a,b,c,f);    
        end
        t = t+dt/2;
        
        u1 = (-cos(pi*x(:))*sin(pi*y(:))')*exp(-2*pi*pi*t/Re);
        u2 = (sin(pi*x(:))*cos(pi*y(:))')*exp(-2*pi*pi*t/Re);
        %u1 = zeros(numNos,numNos);                                           %MUDAR
        %u2 = zeros(numNos,numNos);        
        rho1 = u1*dt/h;
        rho2 = u2*dt/h;
        
        % condicoes de contorno
        for cont = 1:numNos
            %X:
            w(1,cont) = 2*pi*cos(pi*x(1))*cos(pi*y(cont))*exp(-2*pi*pi*t/Re);
            w(numNos,cont) = 2*pi*cos(pi*x(numNos))*cos(pi*y(cont))*exp(-2*pi*pi*t/Re);                
            %Y:
            w(cont,1) = 2*pi*cos(pi*x(cont))*cos(pi*y(1))*exp(-2*pi*pi*t/Re);
            w(cont,numNos) = 2*pi*cos(pi*x(cont))*cos(pi*y(numNos))*exp(-2*pi*pi*t/Re);
        end
        
        %montando o termo fonte para o segundo passo do ADI implicito em Y
        for i = 2:numNos-1
            for j = 2:numNos-1
                a(1,j-1) = -1*(sigma/2+rho2(i,j)/4);%diagonal inferior
                b(1,j-1) = (1+sigma); %diagonal principal
                c(1,j-1) = (-sigma/2+rho2(i,j)/4);%diagonal superior    
                %notar que o rho1 e rho2 variam de 2 até numNos-1 apenas
                f(j-1)= (1-sigma)*w_meio(i,j)+ (sigma/2-(rho1(i,j)/4))*w_meio(i+1,j)+ (sigma/2+(rho1(i,j)/4))*w_meio(i-1,j);
            end
            %COLOCAR AQUI OS TERMOS FONTE QDO i=1 E i=4
            f(1) = f(1)+ w(i,1)*(sigma/2+rho2(i,2)/4); %rho2(i,2) faz parte do contorno
            f(numNos-2) = f(numNos-2) - (-sigma/2+rho2(i,numNos-1)/4)*w(i,numNos); %rho2(i,numNos-1) faz parte do contorno
            w(i,2:numNos-1) = TDMAsolver(a,b,c,f);    
        end
    end
    t = tempoFinal;
    for i = 1:numNos
        for j = 1:numNos
            %matErro(i,j) = abs(w(i,j)-exp(-pi*pi*t*2)*sin(pi*x(i))*sin(pi*y(j)));
            matErro(i,j) = abs(w(i,j)-(2*pi*cos(pi*x(i))*cos(pi*y(j))*exp(-2*pi*pi*t/Re)));     %MUDAR
        end
    end
    erro(k) = max(max(matErro));
    dh(k) = h;
end

plot(-log(dh), log(erro));
xlabel("-log(h)");
ylabel("log(erro)");

figure
r = 2*pi*cos(pi*x(:))*cos(pi*y(:))'*exp(-2*pi*pi*0/Re);
surf(x,y,w);
xlabel("X");
ylabel("Y");
zlabel("\omega");
title("Solução aproximada");
figure
%surf(x,y,exp(-pi*pi*t*2)*sin(pi*x(:))*sin(pi*y(:)'));                     %MUDAR
surf(x,y,2*pi*cos(pi*x(:))*cos(pi*y(:))'*exp(-2*pi*pi*t/Re));
xlabel("X");
ylabel("Y");
zlabel("\omega");
title("Solução exata");

% figure
% surf(x,y,r);
% xlabel("X");
% ylabel("Y");
% zlabel("\omega");
% title("Condição inicial")