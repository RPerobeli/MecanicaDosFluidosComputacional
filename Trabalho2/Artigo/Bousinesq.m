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
numElementos = 31;
numNos = numElementos +1;
%Dimensão do elemento
h = (Xf-Xi)/numElementos;
h2 = h*h;
%passo temporal
dt = h2;
%numero de Reynolds, Rayleigh e Prandtl
Re = 1000;
Ra = 1.e4;
Pr = 0.7;
%alfas
alfa_psi=1;
alfa_t=1;
alfa_w=1;
%coordenadas:
x= Xi:h:Xf;
y= Yi:h:Xf;

%Condição inicial
w = zeros(numNos,numNos);
w_meio = zeros(numNos,numNos);
psi = zeros(numNos,numNos);
psi_meio = zeros(numNos,numNos);
T = zeros(numNos,numNos);
T_meio = zeros(numNos,numNos);

sigma = dt/(2*h2);

rhox = zeros(numNos,numNos);
rhoy = zeros(numNos, numNos);
theta_x = zeros(numNos, numNos);

Residuos = ones(3,1);
critParada = 1;

while (critParada > 1e-1)
    %% Função Corrente    

    %METODO ADI
    %primeiro passo do ADI implicito em X
    for j = 2:numNos-1
        a = -sigma*ones(1,numNos-2);%diagonal inferior
        b = ((1./alfa_psi)+2.*sigma)*ones(1,numNos-2); %diagonal principal
        c = -sigma*ones(1,numNos-2);%diagonal superior 
        for i = 2:numNos-1  
            f(i-1)= w(i,j)*(dt/2) + (sigma)*psi(i,j+1)+((1/alfa_psi)-2*sigma)*psi(i,j) +(sigma)*psi(i,j-1);
        end
        %COLOCAR AQUI OS TERMOS FONTE QDO i=1 E i=4
%         f(1) = f(1)+ psi_meio(1,j)*sigma;
%         f(numNos-2) = f(numNos-2) + (sigma)*psi_meio(numNos,j);
        psi_meio(2:numNos-1,j) = TDMAsolver(a,b,c,f);    
    end

    %montando o termo fonte para o segundo passo do ADI implicito em Y
    for i = 2:numNos-1
        a = -sigma*ones(1,numNos-2);%diagonal inferior
        b = ((1./alfa_psi)+2*sigma)*ones(1,numNos-2); %diagonal principal
        c = -sigma*ones(1,numNos-2);%diagonal superior 
        for j = 2:numNos-1
            f(j-1)= w(i,j)*(dt/2) + (sigma)*psi_meio(i+1,j)+((1./alfa_psi)-2*sigma)*psi_meio(i,j)+(sigma)*psi_meio(i-1,j);
        end
        %COLOCAR AQUI OS TERMOS FONTE QDO i=1 E i=4
%         f(1) = f(1)+ psi(i,1)*(sigma);
%         f(numNos-2) = f(numNos-2)+(sigma)*psi(i,numNos);
        psi(i,2:numNos-1) = TDMAsolver(a,b,c,f);    
    end
    
    %% Temperatura
    %Condições de Contorno (Neumann em Y)
    %Y:
    %T_meio(:,1) = 0;
    %T_meio(:,numNos) = 0; 
    %X:
    T_meio(numNos,:) = 0;
    T_meio(1,:) = 1; %esquerda por ultimo
    
    for i = 2:numNos-1
            rhox(i,:) = 1/4 * (psi(i+1,:)-psi(i-1,:));
            rhoy(:,i) = 1/4 * (psi(:,i+1) - psi(:,i-1));
    end
    
    %Metodo ADI
    %implicito em X:
    for j= 1:numNos
        if(j==1)
            for i = 2:numNos-1
                %notar que o rhox e rhoy variam de 2 até numNos-1            
                a(1,i-1) = -sigma*(rhoy(i,j+1)+1);%diagonal inferior
                b(1,i-1) = (1./alfa_t)+2.*sigma; %diagonal principal
                c(1,i-1) = sigma*(rhoy(i,j+1)-1); %diagonal superior

                f(i-1) = ((1./alfa_t)-2.*sigma)*T(i,j)+sigma*(rhox(i,j+1)+1)*T(i,j+1) + sigma*(1-rhox(i,j+1))*T(i,j+1);
            end
        else if(j == numNos)
                for i = 2:numNos-1
                    %notar que o rhox e rhoy variam de 2 até numNos-1            
                    a(1,i-1) = -sigma*(rhoy(i,j-1)+1);%diagonal inferior
                    b(1,i-1) = (1./alfa_t)+2.*sigma; %diagonal principal
                    c(1,i-1) = sigma*(rhoy(i,j-1)-1); %diagonal superior

                    f(i-1) = ((1./alfa_t)-2*sigma)*T(i,j)+sigma*(rhox(i,j-1)+1)*T(i,j-1) + sigma*(1-rhox(i,j-1))*T(i,j-1);
                end
            else
                for i = 2:numNos-1
                    %notar que o rhox e rhoy variam de 2 até numNos-1            
                    a(1,i-1) = -sigma*(rhoy(i,j)+1);%diagonal inferior
                    b(1,i-1) = (1./alfa_t)+2*sigma; %diagonal principal
                    c(1,i-1) = sigma*(rhoy(i,j)-1); %diagonal superior

                    f(i-1) = ((1./alfa_t)-2*sigma)*T(i,j)+sigma*(rhox(i,j)+1)*T(i,j+1) + sigma*(1-rhox(i,j))*T(i,j-1);
                end
            end
            
        end
        %COLOCAR AQUI OS TERMOS FONTE QDO i=1 E i=4
        f(1) = f(1)+ T_meio(1,j)*sigma*(rhoy(2,j)+1);
        f(numNos-2) = f(numNos-2) - sigma*(rhoy(numNos-1,j)-1)*T_meio(numNos,j);
        T_meio(2:numNos-1,j) = TDMAsolver(a,b,c,f);
    end
    
    %condiçoes de contorno: por mais que não variem no tempo, é necessário
    %para colocar para alterar a matriz de vorticidade final
    %Y:
%     T(:,1) = 0;
%     T(:,numNos) = 0; 
    %X:
    T(numNos,:) = 0;
    T(1,:) = 1;%esquerda por ultimo
    
    %Implicito em Y:
   % dT/dy = 0 => (T_{i,j+1} - T_{i,j-1})/2h = 0 => T_{i,j+1} = T_{i,j-1}, 
    for i= 2:numNos-1
        for j = 1:numNos
            if(j ==1)
                %notar que o rhox e rhoy variam de 2 até numNos-1
                o(1,j) = 0;%diagonal inferior
                p(1,j) = (1./alfa_t)+2*sigma; %diagonal principal
                q(1,j) = -sigma*2; %diagonal superior
                r(j) = ((1./alfa_t)-2*sigma)*T_meio(i,j)+sigma*(1-rhoy(i,j+1))*T_meio(i+1,j) + sigma*(1+rhoy(i,j+1))*T_meio(i-1,j);
            else if (j == numNos)
                 %notar que o rhox e rhoy variam de 2 até numNos-1
                 o(1,j) = -sigma*2;%diagonal inferior
                 p(1,j) = (1./alfa_t)+2*sigma; %diagonal principal
                 q(1,j) = 0; %diagonal superior
                 r(j) = ((1./alfa_t)-2*sigma)*T_meio(i,j)+sigma*(1-rhoy(i,j-1))*T_meio(i+1,j) + sigma*(1+rhoy(i,j-1))*T_meio(i-1,j);
                 else
                    %notar que o rhox e rhoy variam de 2 até numNos-1
                    o(1,j) = sigma*(rhox(i,j)-1);%diagonal inferior
                    p(1,j) = (1./alfa_t)+2*sigma; %diagonal principal
                    q(1,j) = -sigma*(rhox(i,j)+1); %diagonal superior 
                    r(j) = ((1./alfa_t)-2*sigma)*T_meio(i,j)+sigma*(1-rhoy(i,j))*T_meio(i+1,j) + sigma*(1+rhoy(i,j))*T_meio(i-1,j);
                 end
            end
%             r(j) = (1-2*sigma)*T_meio(i,j)+sigma*(1-rhoy(i,j))*T_meio(i+1,j) + sigma*(1+rhoy(i,j))*T_meio(i-1,j);
        end
        %COLOCAR AQUI OS TERMOS FONTE QDO i=1 E i=4
        %r(1) = r(1) - T(i,1)*sigma*(rhox(i,2)-1);
        %r(numNos) = r(numNos) + sigma*(rhox(i,numNos)+1)*T(i,numNos);
        T(i,1:numNos) = TDMAsolver(o,p,q,r);
    end

    %% Vorticidade
    %Condições de contorno
    %X:
    w_meio(1,:) = -2*psi(2,:)/h2;
    w_meio(numNos,:) = -2*psi(numNos-1,:)/h2;                
    %Y:
    w_meio(:,1) = -2*psi(:,2)/h2;
    w_meio(:,numNos) = -2*psi(:,numNos-1)/h2; %- 2/h; %superior por ultimo
    
    for i = 2:numNos-1
            rhox(i,:) = 1/4 * (psi(i+1,:)-psi(i-1,:));
            rhoy(:,i) = 1/4 * (psi(:,i+1) - psi(:,i-1));
            theta_x(i,:)= (alfa_w*dt)/(4*h)*(T(i+1,:)-T(i-1,:));
    end
    
    %METODO ADI
    %Implícito em X:    
    for j= 2:numNos-1
%             a = -sigma*(rhoy(2:numNos-1,j)+Pr);%diagonal inferior
%             b = (1+2*Pr*sigma)*ones(1,numNos-2); %diagonal principal
%             c = sigma*(rhoy(2:numNos-1,j)-Pr);%diagonal superior 

        for i = 2:numNos-1
            %notar que o rhox e rhoy variam de 2 até numNos-1            
            a(1,i-1) = -sigma*(rhoy(i,j)+Pr);%diagonal inferior
            b(1,i-1) = (1./alfa_w)+2*Pr*sigma; %diagonal principal
            c(1,i-1) = sigma*(rhoy(i,j)-Pr); %diagonal superior
            
            f(i-1) = ((1./alfa_w)-2*Pr*sigma)*w(i,j)+sigma*(rhox(i,j)+Pr)*w(i,j+1) + sigma*(Pr-rhox(i,j))*w(i,j-1)+ Ra*Pr*theta_x(i,j);
        end
        %COLOCAR AQUI OS TERMOS FONTE QDO i=1 E i=4
        f(1) = f(1)+ w_meio(1,j)*sigma*(rhoy(2,j)+Pr);
        f(numNos-2) = f(numNos-2) - sigma*(rhoy(numNos-1,j)-Pr)*w_meio(numNos,j);
        w_meio(2:numNos-1,j) = TDMAsolver(a,b,c,f);
    end

    %condiçoes de contorno: por mais que não variem no tempo, é necessário
    %para colocar para alterar a matriz de vorticidade final
    %X:
    w(1,:) = -2*psi(2,:)/h2;
    w(numNos,:) = -2*psi(numNos-1,:)/h2;                
    %Y:
    w(:,1) = -2*psi(:,2)/h2;
    w(:,numNos) = -2*psi(:,numNos-1)/h2;% - 2/h; %superior por ultimo
    
    %Implicito em Y:
    
    for i= 2:numNos-1
        for j = 2:numNos-1
            %notar que o rhox e rhoy variam de 2 até numNos-1
            a(1,j-1) = sigma*(rhox(i,j)-Pr);%diagonal inferior
            b(1,j-1) = (1./alfa_w)+2*Pr*sigma; %diagonal principal
            c(1,j-1) = -sigma*(rhox(i,j)+Pr); %diagonal superior
            
            f(j-1) = ((1./alfa_w)-2*Pr*sigma)*w_meio(i,j)+sigma*(Pr-rhoy(i,j))*w_meio(i+1,j) + sigma*(Pr+rhoy(i,j))*w_meio(i-1,j) + Ra*Pr*theta_x(i,j);
        end
        %COLOCAR AQUI OS TERMOS FONTE QDO i=1 E i=4
        f(1) = f(1) - w(i,1)*sigma*(rhox(i,2)-Pr);
        f(numNos-2) = f(numNos-2) + sigma*(rhox(i,numNos-1)+Pr)*w(i,numNos);
        w(i,2:numNos-1) = TDMAsolver(a,b,c,f);
    end
    
    for i = 2:numNos-1
        for j= 2:numNos-1
            Residuo_psi(i-1,j-1)  = (psi(i+1,j)-2*psi(i,j)+psi(i-1,j))/h2 + (psi(i,j+1)-2*psi(i,j)+psi(i,j-1))/h2 + w(i,j);
            
            A = ((psi(i,j+1)-psi(i,j-1))/(2*h))*((w(i+1,j)-w(i-1,j))/(2*h));
            B = ((psi(i+1,j)-psi(i-1,j))/(2*h))*((w(i,j+1)-w(i,j-1))/(2*h));
            A1 = Ra*Pr*(T(i+1,j)-T(i-1,j))/(2*h);
            Residuo_w(i-1,j-1) = Pr*((w(i+1,j)-2*w(i,j)+w(i-1,j))/(h2) + (w(i,j+1)-2*w(i,j)+w(i,j-1))/(h2)) - A + B + A1;
            
            C = ((psi(i,j+1)-psi(i,j-1))/(2*h))*((T(i+1,j)-T(i-1,j))/(2*h));
            D = ((psi(i+1,j)-psi(i-1,j))/(2*h))*((T(i,j+1)-T(i,j-1))/(2*h));
            Residuo_T(i-1,j-1) = (T(i+1,j)-2*T(i,j)+T(i-1,j))/(h2) + (T(i,j+1)-2*T(i,j)+T(i,j-1))/(h2) - C + D;
        end
    end
    %Residuos(1) = Residuos(1) - 0.25;
    %Residuos(2) = 0;
    Residuos(1)= max(max(abs(Residuo_psi)));
    Residuos(2)=max(max(abs(Residuo_w)));
    Residuos(3)= max(max(abs(Residuo_T)));
    critParada = max(Residuos)
end

for i=1:numNos 
    dT(i) = (T(3,i) - T(1,i))/(2.*h); 
end
Nu0 = 0.;
for i=1:numNos-1 %integração 
    Nu0 = Nu0 + (dT(i)+dT(i+1))*h/2.;
end

maxPsi = max(max(abs(psi)));
maxW = max(max(abs(w)));
figure
[C,h]=contour(x,y,psi','k');
clabel(C,h);
xlabel("X");
ylabel("Y");
zlabel("\Psi");
title("Função Corrente Pr="+Pr+ "Ra ="+Ra);%+ "\Psi ="+maxPsi);

figure
[C,h]=contour(x,y,w','k');
clabel(C,h);
xlabel("X");
ylabel("Y");
zlabel("\omega");
title("Vorticidade Pr="+Pr+ "Ra ="+Ra);

figure
[C,h]=contour(x,y,T','k');
clabel(C,h);
xlabel("X");
ylabel("Y");
zlabel("T");
title("Temperatura Pr="+Pr+ "Ra ="+Ra);

