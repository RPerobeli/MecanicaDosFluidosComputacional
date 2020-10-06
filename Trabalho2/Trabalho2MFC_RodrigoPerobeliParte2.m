clear all
close all
clc

%Dominio temporal:
tempoInicial = 0.0;             
%Dominio espacial
Xi = 0.0;
Xf = 1.0;                      
Yi = 0.0;
Yf = 1.0;
numElementos = 32;
numNos = numElementos +1;
%Dimensão do elemento
h = (Xf-Xi)/numElementos;
h2 = h*h;
%passo temporal
dt = h2/10;
%numero de Reynolds
Re = 1000;
%coordenadas:
x= Xi:h:Xf;
y= Yi:h:Xf;
% x = zeros(numNos,1);
% y = zeros(numNos,1);
% x(1) = Xi;
% y(1) = Yi;
% for i = 2:numNos
%     %Define um vetor de coordenadas 
%     x(i) = Xi+(i-1)*h;
%     y(i) = Yi+(i-1)*h;
% end
%Condição inicial
t= tempoInicial;
w = zeros(numNos,numNos);
w_meio = zeros(numNos,numNos);
psi = zeros(numNos,numNos);
psi_meio = zeros(numNos,numNos);
sigma = dt/(2*h2);

rhox = zeros(numNos,numNos);
rhoy = zeros(numNos, numNos);

Residuos = ones(2,1);
critParada = 1;

while (critParada > 1e-6)
    %% Função Corrente
    %condições de contorno
    %X:
%     psi_meio(1,:) = 0;
%     psi_meio(numNos,:) = 0;                
%     %Y:
%     psi_meio(:,1) = 0;
%     psi_meio(:,numNos) = 0;
    
    %METODO ADI
    %primeiro passo do ADI implicito em X
    for j = 2:numNos-1
        a = -sigma*ones(1,numNos-2);%diagonal inferior
        b = (1+2*sigma)*ones(1,numNos-2); %diagonal principal
        c = -sigma*ones(1,numNos-2);%diagonal superior 
        for i = 2:numNos-1  
            f(i-1)= w(i,j)*(dt/2) + (sigma)*psi(i,j+1)+(1-2*sigma)*psi(i,j) +(sigma)*psi(i,j-1);
        end
        %COLOCAR AQUI OS TERMOS FONTE QDO i=1 E i=4
%         f(1) = f(1)+ psi_meio(1,j)*sigma;
%         f(numNos-2) = f(numNos-2) + (sigma)*psi_meio(numNos,j);
        psi_meio(2:numNos-1,j) = TDMAsolver(a,b,c,f);    
    end

    % condicoes de contorno
        %X:
%         psi(1,:) = 0;
%         psi(numNos,:) = 0;                
%         %Y:
%         psi(:,1) = 0;
%         psi(:,numNos) = 0;
    
    %montando o termo fonte para o segundo passo do ADI implicito em Y
    for i = 2:numNos-1
        a = -sigma*ones(1,numNos-2);%diagonal inferior
        b = (1+2*sigma)*ones(1,numNos-2); %diagonal principal
        c = -sigma*ones(1,numNos-2);%diagonal superior 
        for j = 2:numNos-1
            f(j-1)= w(i,j)*(dt/2) + (sigma)*psi_meio(i+1,j)+(1-2*sigma)*psi_meio(i,j)+(sigma)*psi_meio(i-1,j);
        end
        %COLOCAR AQUI OS TERMOS FONTE QDO i=1 E i=4
%         f(1) = f(1)+ psi(i,1)*(sigma);
%         f(numNos-2) = f(numNos-2)+(sigma)*psi(i,numNos);
        psi(i,2:numNos-1) = TDMAsolver(a,b,c,f);    
    end
    
    %% Vorticidade    
    %Condições de contorno
    %X:
    w_meio(1,:) = -2*psi(2,:)/h2;
    w_meio(numNos,:) = -2*psi(numNos-1,:)/h2;                
    %Y:
    w_meio(:,1) = -2*psi(:,2)/h2;
    w_meio(:,numNos) = -2*psi(:,numNos-1)/h2 - 2/h; %superior por ultimo
    
    for i = 2:numNos-1
            rhox(i,:) = Re/4 * (psi(i+1,:)-psi(i-1,:));
            rhoy(:,i) = Re/4 * (psi(:,i+1) - psi(:,i-1));
    end
    
    %METODO ADI
    %Implícito em X:    
    for j= 2:numNos-1
        for i = 2:numNos-1
            %notar que o rhox e rhoy variam de 2 até numNos-1            
            a(1,i-1) = -sigma*(rhoy(i,j)+1);%diagonal inferior
            b(1,i-1) = 1+2*sigma; %diagonal principal
            c(1,i-1) = sigma*(rhoy(i,j)-1); %diagonal superior
            
            f(i-1) = (1-2*sigma)*w(i,j)+sigma*(rhox(i,j)+1)*w(i,j+1) + sigma*(1-rhox(i,j))*w(i,j-1);
        end
        %COLOCAR AQUI OS TERMOS FONTE QDO i=1 E i=4
        f(1) = f(1)+ w_meio(1,j)*sigma*(rhoy(2,j)+1);
        f(numNos-2) = f(numNos-2) - sigma*(rhoy(numNos-1,j)-1)*w_meio(numNos,j);
        w_meio(2:numNos-1,j) = TDMAsolver(a,b,c,f);
    end
    
    %condiçoes de contorno: por mais que não variem no tempo, é necessário
    %para colocar para alterar a matriz de vorticidade final
    %X:
    w(1,:) = -2*psi(2,:)/h2;
    w(numNos,:) = -2*psi(numNos-1,:)/h2;                
    %Y:
    w(:,1) = -2*psi(:,2)/h2;
    w(:,numNos) = -2*psi(:,numNos-1)/h2 - 2/h; %superior por ultimo
    
    %Implicito em Y:
    
    for i= 2:numNos-1
        for j = 2:numNos-1
            %notar que o rhox e rhoy variam de 2 até numNos-1
            a(1,j-1) = sigma*(rhox(i,j)-1);%diagonal inferior
            b(1,j-1) = 1+2*sigma; %diagonal principal
            c(1,j-1) = -sigma*(rhox(i,j)+1); %diagonal superior
            
            f(j-1) = (1-2*sigma)*w_meio(i,j)+sigma*(1-rhoy(i,j))*w_meio(i+1,j) + sigma*(1+rhoy(i,j))*w_meio(i-1,j);
        end
        %COLOCAR AQUI OS TERMOS FONTE QDO i=1 E i=4
        f(1) = f(1) - w(i,1)*sigma*(rhox(i,2)-1);
        f(numNos-2) = f(numNos-2) + sigma*(rhox(i,numNos-1)+1)*w(i,numNos);
        w(i,2:numNos-1) = TDMAsolver(a,b,c,f);
    end
    
    for i = 2:numNos-1
        for j= 2:numNos-1
            Residuo_psi(i-1,j-1)  = (psi(i+1,j)-2*psi(i,j)+psi(i-1,j))/h2 + (psi(i,j+1)-2*psi(i,j)+psi(i,j-1))/h2 + w(i,j);
            A = ((psi(i,j+1)-psi(i,j-1))/(2*h))*((w(i+1,j)-w(i-1,j))/(2*h));
            B = ((psi(i+1,j)-psi(i-1,j))/(2*h))*((w(i,j+1)-w(i,j-1))/(2*h));
            Residuo_w(i-1,j-1) = (w(i+1,j)-2*w(i,j)+w(i-1,j))/(Re*h2) + (w(i,j+1)-2*w(i,j)+w(i,j-1))/(Re*h2) - A + B;
        end
    end
    %Residuos(1) = Residuos(1) - 0.25;
    %Residuos(2) = 0;
    Residuos(1)= max(max(abs(Residuo_psi)));
    Residuos(2)=max(max(abs(Residuo_w)));
    critParada = max(Residuos)
end

maxPsi = max(max(abs(psi)));
maxW = max(max(abs(w)));
figure
surf(x,y,psi);
xlabel("X");
ylabel("Y");
zlabel("\Psi");
title("Função Corrente Re="+Re+ "\Psi ="+maxPsi);
figure
surf(x,y,w);
xlabel("X");
ylabel("Y");
zlabel("\omega");
title("Vorticidade Re="+Re);


