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

%for k =1:NumMalhas
    %O número de elementos deve ser definido de forma que possa ser variável
    %para depois fazer a convergência
    %numElementos = 2^(k+1);
    numElementos = 32;
    numNos = numElementos+1;
    %Dimensão do elemento
    h = (Xf-Xi)/numElementos;
    h2 = h*h;
    %passo temporal
    dt = h;
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
            %Matriz que será solução intermediaria
            w_meio(i,j)=w(i,j);
        end
    end
    %sigma é um numero de CFL,
    Pe=1/Re; %Pe = epsilon = numero de pecle = 1
    sigma = Pe*(dt/h2);
    
    while t < tempoFinal
        %resolve para cada tempo a malha inteira
        t = t+dt/2;
        
        %velocidade de convecção
        u1 = -cos(pi*x(:))*sin(pi*y(:)')*exp(-2*pi*t/Re);                  %MUDAR
        u2 = sin(pi*x(:))*cos(pi*y(:)')*exp(-2*pi*t/Re);
        %u1 = zeros(numNos,1);
        %u2 = zeros(numNos,1);
        
        rho1 = u1*dt/h;
        rho2 = u2*dt/h;
        
        %METODO ADI
        %primeiro passo do ADI implicito em X
        for j = 2:numNos-1
            a = -1*(sigma/2+(rho1(j)/4))*ones(1, numNos-2);%diagonal inferior
            b = (1+sigma)*ones(1, numNos-2); %diagonal principal
            c = ((-sigma/2)+(rho1(j)/4))*ones(1, numNos-2);%diagonal superior           
            for i = 2:numNos-1
                f(i-1)= (1-sigma)*w(i,j)+ (sigma/2-rho2(j)/4)*w(i,j+1)+ ((sigma/2)+(rho2(j)/4))*w(i,j-1);
            end
            %COLOCAR AQUI OS TERMOS FONTE QDO i=1 E i=4
            %w(1,j) = exp(-pi*pi*t*2)*sin(pi*x(1))*sin(pi*y(j));
            %w(numNos,j) = exp(-pi*pi*t*2)*sin(pi*x(numNos))*sin(pi*y(j));
            f(1) = f(1)+ w(1,j)*((sigma/2)+(rho1(j)/4));
            f(numNos-2) = f(numNos-2) - ((-sigma/2)+(rho1(j)/4))*w(numNos,j);
            w_meio(2:numNos-1,j) = TDMAsolver(a,b,c,f);    
        end
        t = t+dt/2;
        
        u1 = -cos(pi*x(:))*sin(pi*y(:)')*exp(-2*pi*t/Re);
        u2 = sin(pi*x(:))*cos(pi*y(:)')*exp(-2*pi*t/Re);
        %u1 = zeros(numNos,1);                                           %MUDAR
        %u2 = zeros(numNos,1);
        
        rho1 = u1*dt/h;
        rho2 = u2*dt/h;
        
        %montando o termo fonte para o segundo passo do ADI implicito em Y
        for i = 2:numNos-1
            a = -1*(sigma/2+rho2(i)/4)*ones(1, numNos-2);%diagonal inferior
            b = (1+sigma)*ones(1, numNos-2); %diagonal principal
            c = (-sigma/2+rho2(i)/4)*ones(1, numNos-2);%diagonal superior    
            for j = 2:numNos-1
                f(j-1)= (1-sigma)*w_meio(i,j)+ (sigma/2-(rho1(i)/4))*w_meio(i+1,j)+ (sigma/2+(rho1(i)/4))*w_meio(i-1,j);
            end
            %COLOCAR AQUI OS TERMOS FONTE QDO i=1 E i=4
            %w_meio(i,1) = exp(-pi*pi*t*2)*sin(pi*x(i))*sin(pi*y(1));
            %w_meio(i,numNos) = exp(-pi*pi*t*2)*sin(pi*x(i))*sin(pi*y(numNos));
            
            f(1) = f(1)+ w_meio(i,1)*(sigma/2+rho2(i)/4);
            f(numNos-2) = f(numNos-2) - (-sigma/2+rho2(i)/4)*w_meio(i,numNos);
            w(i,2:numNos-1) = TDMAsolver(a,b,c,f);    
        end
    end
    t = tempoFinal;
    for i = 1:numNos
        for j = 1:numNos
            %matErro(i,j) = abs(w(i,j)-exp(-pi*pi*t*2)*sin(pi*x(i))*sin(pi*y(j)));
            matErro(i,j) = abs(w(i,j)-2*pi*cos(pi*x(i))*cos(pi*y(j))*exp(-2*pi*t/Re));     %MUDAR
        end
    end
    %erro(k) = max(max(matErro));
    %dh(k) = h;
%end

%plot(-log(dh), log(erro));
figure
surf(x,y,w);
figure
%surf(x,y,exp(-pi*pi*t*2)*sin(pi*x(:))*sin(pi*y(:)'));                     %MUDAR
surf(x,y,2*pi*cos(pi*x(:))*cos(pi*y(:)')*exp(-2*pi*t/Re));
  
  