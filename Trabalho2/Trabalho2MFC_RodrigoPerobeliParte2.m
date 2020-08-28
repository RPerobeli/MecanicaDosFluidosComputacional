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
numElementos = 4;
numNos = numElementos +1;
%Dimensão do elemento
h = (Xf-Xi)/numElementos;
h2 = h*h;
%passo temporal
dt = h;
%numero de Reynolds
Re = 20;
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
w = zeros(numNos,numNos);
psi = zeros(numNos,numNos);


