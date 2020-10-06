% programa diferencas finitas 
clear all
close all
clc
% domínio espacial [xi,xf]x[xi,xf]
xi = 0.;
xf = 1.;

% Prandtl
Pr = 0.7;

% Rayleigh
Ra = 1.e5;

% numero de elementos 
nel = 81; % numero de nos = nel + 1

% dimensao do elemento
h = (xf-xi)/nel;
h2 = h*h;
dt = h2;

% coordenadas
x = xi:h:xf;
y = xi:h:xf;

sig = 0.5*dt/h2;

alpha_T   = 0.5;
alpha_psi = 0.1;
alpha_w   = 0.5;

% definição dos vetores de variaveis
% stream function
psix = zeros(nel+1,nel+1);
psiy = zeros(nel+1,nel+1);
% vorticity
wx = zeros(nel+1,nel+1);
wy = zeros(nel+1,nel+1);
% temperature
Tx = zeros(nel+1,nel+1);
Ty = zeros(nel+1,nel+1); 

t = 0;
tf = 1000;
error = 1.;

while error > 1.e-1
    %% step 1 stream function
    for j=2:nel
        % montagem da matriz
        a = -sig*ones(1,nel-1); % diagonal inferior
        b = (1./alpha_psi + 2.*sig)*ones(1,nel-1);  % diagonal principal
        c = -sig*ones(1,nel-1); % diagonal superior

        % montagem do vetor fonte
        f = dt*wy(2:nel,j)/2. + (1./alpha_psi-2.*sig)*psiy(2:nel,j) ...
               + sig*(psiy(2:nel,j+1) + psiy(2:nel,j-1));
           
        % resolucao do sistema pelo algoritmo de Thomas
        psix(2:nel,j) = TDMAsolver(a,b,c,f);
    end   
    %step 2 stream function
    for i=2:nel
        % montagem da matriz
        a = -sig*ones(1,nel-1); % diagonal inferior
        b = (1./alpha_psi + 2.*sig)*ones(1,nel-1);  % diagonal principal
        c = -sig*ones(1,nel-1); % diagonal superior

        % montagem do vetor fonte
        f = dt*wy(i,2:nel)/2. + (1./alpha_psi-2.*sig)*psix(i,2:nel) ...
               + sig*(psix(i+1,2:nel) + psix(i-1,2:nel));
           
        % resolucao do sistema pelo algoritmo de Thomas
        psiy(i,2:nel) = TDMAsolver(a,b,c,f);
    end 
    
    %% Temp
    % condicoes de contorno
    Tx(1,:) = 1.; % left
    Tx(nel+1,:) = 0.; % right
    Ty(1,:) = 1.; % left
    Ty(nel+1,:) = 0.; % right


    for i=2:nel
        rhox(i,:) = (psiy(i+1,:) - psiy(i-1,:))/(4.);
        rhoy(:,i) = (psiy(:,i+1) - psiy(:,i-1))/(4.);
    end
        rhox(1,:) = rhox(2,:);
        rhoy(:,1) = rhoy(:,2);
        rhox(nel+1,:) = rhox(nel,:);
        rhoy(:,nel+1) = rhoy(:,nel);
        
        
    %step 1 temperature
    for j=1:nel+1
        % montagem da matriz
        a = -sig*(rhoy(2:nel,j) + 1.); % diagonal inferior
        b = (1./alpha_T + 2.*sig)*ones(1,nel-1);  % diagonal principal
        c =  sig*(rhoy(2:nel,j) - 1.); % diagonal superior

        % montagem do vetor fonte
        for i = 2:nel
            if(j==1) 
                f(i-1) = (1./alpha_T-2.*sig)*Ty(i,j) + sig*(rhox(i,j)+1.)*Ty(i,j+1) ...
                                             + sig*(1.-rhox(i,j))*Ty(i,j+1);
            else if(j==nel+1) 
                f(i-1) = (1./alpha_T-2.*sig)*Ty(i,j) + sig*(rhox(i,j)+1.)*Ty(i,j-1) ...
                                             + sig*(1.-rhox(i,j))*Ty(i,j-1); 
                else
                f(i-1) = (1./alpha_T-2.*sig)*Ty(i,j) + sig*(rhox(i,j)+1.)*Ty(i,j+1) ...
                                             + sig*(1.-rhox(i,j))*Ty(i,j-1); 
                end
            end
        end 
    % tratamento condicoes de contorno
        f(1) = f(1) + sig*(rhoy(2,j) + 1.)*Tx(1,j);
        f(nel-1) = f(nel-1) - sig*(rhoy(nel,j) - 1.)*Tx(nel+1,j);
        % resolucao do sistema pelo algoritmo de Thomas
        Tx(2:nel,j) = TDMAsolver(a,b,c,f);
    end   
    %step 2 temperature
    for i=2:nel
         % montagem da matriz
        a1 = sig*(rhox(i,1:nel+1) - 1.); % diagonal inferior
        b1 = (1./alpha_T + 2.*sig)*ones(1,nel+1);  % diagonal principal
        c1 = -sig*(rhox(i,1:nel+1) + 1.); % diagonal superior
        
        c1(1) = c1(1) + sig*(rhox(i,1) - 1.);
        a1(nel+1) = a1(nel+1) - sig*(rhox(i,nel+1) + 1.);
        
        % montagem do vetor fonte
        for j = 1:nel+1
            f1(j) = (1./alpha_T-2.*sig)*Tx(i,j) + sig*(1.-rhoy(i,j))*Tx(i+1,j) ...
                                       + sig*(1.+rhoy(i,j))*Tx(i-1,j);
        end 
           
        % resolucao do sistema pelo algoritmo de Thomas
        Ty(i,1:nel+1) = TDMAsolver(a1,b1,c1,f1);
    end 

    %% Vort
   % condicoes de contorno
    wx(1,:) = -2.*psiy(2,:)/h2; % left
    wx(nel+1,:) = -2.*psiy(nel,:)/h2; % right
    wx(:,1) = -2.*psiy(:,2)/h2; % bottom
    wx(:,nel+1) =  -2.*psiy(:,nel)/h2; % top
    
    for i=2:nel
        theta(i,:) = dt*Ra*Pr*(Ty(i+1,:)-Ty(i-1,:))/(4.*h);
    end
    % primeiro passo solução em x
    for j=2:nel
        % montagem da matriz
        a = -sig*(rhoy(2:nel,j) + Pr); % diagonal inferior
        b = (1./alpha_w + 2.*Pr*sig)*ones(1,nel-1);  % diagonal principal
        c =  sig*(rhoy(2:nel,j) - Pr); % diagonal superior

        % montagem do vetor fonte
        for i = 2:nel
            f(i-1) = theta(i,j) + (1./alpha_w-2.*Pr*sig)*wy(i,j) + sig*(rhox(i,j)+Pr)*wy(i,j+1) ...
                                                         + sig*(Pr-rhox(i,j))*wy(i,j-1);
        end
    % tratamento condicoes de contorno
        f(1) = f(1) + sig*(rhoy(2,j) + Pr)*wx(1,j);
        f(nel-1) = f(nel-1) - sig*(rhoy(nel,j) - Pr)*wx(nel+1,j);

        % resolucao do sistema pelo algoritmo de Thomas
        wx(2:nel,j) = TDMAsolver(a,b,c,f);
    end 

    wy(1,:) = wx(1,:); % left
    wy(nel+1,:) = wx(nel+1,:); % right
    wy(:,1) = wx(:,1); % bottom
    wy(:,nel+1) = wx(:,nel+1); % top

    % segundo passo solução em y       
    for i=2:nel
        % montagem da matriz
        a = sig*(rhox(i,2:nel) - Pr); % diagonal inferior
        b = (1./alpha_w + 2.*Pr*sig)*ones(1,nel-1);  % diagonal principal
        c = -sig*(rhox(i,2:nel) + Pr); % diagonal superior

        % montagem do vetor fonte
        for j = 2:nel
            f(j-1) = theta(i,j) + (1./alpha_w-2.*Pr*sig)*wx(i,j) + sig*(Pr-rhoy(i,j))*wx(i+1,j) ...
                                                        + sig*(Pr+rhoy(i,j))*wx(i-1,j);
        end 

        % tratamento condicoes de contorno
        f(1) = f(1) - sig*(rhox(i,2) - Pr)*wy(i,1);
        f(nel-1) = f(nel-1) + sig*(rhox(i,nel) + Pr)*wy(i,nel+1);

        % resolucao do sistema pelo algoritmo de Thomas
        wy(i,2:nel) = TDMAsolver(a,b,c,f);
    end 
    
    t = t + 1;
    if (mod(t,100) == 0)
        for i=2:nel
            for j=2:nel
                Rpsi(i,j) = (psiy(i+1,j) - 2.*psiy(i,j) + psiy(i-1,j) ...
                          +  psiy(i,j+1) - 2.*psiy(i,j) + psiy(i,j-1))/h2 + wy(i,j);

                Rw(i,j) = Pr*(wy(i+1,j) - 2.*wy(i,j) + wy(i-1,j) ...
                        + wy(i,j+1) - 2.*wy(i,j) + wy(i,j-1))/(Ra*h2) ...
                        - (psiy(i,j+1) - psiy(i,j-1))*(wy(i+1,j) - wy(i-1,j))/(Ra*4.*h2)...
                        + (psiy(i+1,j) - psiy(i-1,j))*(wy(i,j+1) - wy(i,j-1))/(Ra*4.*h2)...
                        + Pr*(Ty(i+1,j)-Ty(i-1,j))/(2.*h);
                RT(i,j) = (Ty(i+1,j) - 2.*Ty(i,j) + Ty(i-1,j) ...
                        + Ty(i,j+1) - 2.*Ty(i,j) + Ty(i,j-1))/h2 ...
                        - (psiy(i,j+1) - psiy(i,j-1))*(Ty(i+1,j) - Ty(i-1,j))/(4.*h2)...
                        + (psiy(i+1,j) - psiy(i-1,j))*(Ty(i,j+1) - Ty(i,j-1))/(4.*h2);
            end
        end
%         error1 = max(max(abs(Rpsi)))
%         error2 = max(max(abs(Rw)))
%         error  = max(error1,error2);
%         error3 = max(max(abs(RT)))
%         error  = max(error,error3);

         error = max(max(max(abs(Rpsi))),max(max(abs(Rw))));
         error = max(error,max(max(abs(RT))))
    end
end

%Nusselt Number
%  for j=2:nel
%     for i=2:nel
%         Q(i-1,j-1) = (psiy(i,j+1) - psiy(i,j-1))*Ty(i,j)/(2.*h) - (Ty(i+1,j) - Ty(i-1,j))/(2.*h);
%     end 
%  end
%  Nu = zeros(nel-1,1);
%  for i=1:nel-1
%     for j=1:nel-2
%         Nu(i) = Nu(i) + (Q(i,j)+Q(i,j+1))*h/2.;
%     end
%  end

for i=1:nel+1 
    dT(i) = (Ty(3,i) - Ty(1,i))/(2.*h); 
end
Nu0 = 0.;
for i=1:nel %integração 
    Nu0 = Nu0 + (dT(i)+dT(i+1))*h/2.;
end

t
% max(max(abs(psiy)))
% figure(1)
% surf(x,y,psiy');
% figure(2)
% clabel(contour(x,y,Ty',[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9]));
% figure(3)
% surf(x,y,wy');  
% figure(4)
% surf(x(2:nel),y(2:nel),Q(1:nel-1,1:nel-1)');  
% figure(5)
% plot(x(2:nel),Nu(1:nel-1));

figure
[C,h]=contour(x,y,psiy','k');
clabel(C,h);
xlabel("X");
ylabel("Y");
zlabel("\Psi");
title("Função Corrente Pr="+Pr+ "Ra ="+Ra);%+ "\Psi ="+maxPsi);

figure
[C,h]=contour(x,y,wy',[-2000, -1000,-700, -500, -200,0, 300, 500,2000],'k');
%[C,h]=contour(x,y,wy',[-700, -500, 0, 300, 500],'k'); 
clabel(C,h);
xlabel("X");
ylabel("Y");
zlabel("\omega");
title("Vorticidade Pr="+Pr+ "Ra ="+Ra);

figure
[C,h]=contour(x,y,wy','k');
clabel(C,h);
xlabel("X");
ylabel("Y");
zlabel("\omega");
title("Vorticidade Pr="+Pr+ "Ra ="+Ra);

figure
[C,h]=contour(x,y,Ty','k');
clabel(C,h);
xlabel("X");
ylabel("Y");
zlabel("T");
title("Temperatura Pr="+Pr+ "Ra ="+Ra);

