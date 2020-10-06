function x = TDMAsolver(a,b,c,d)
% a, b, c são os vetores colunas para a compressão da matriz tridiagonal, d é o vetor à direita
% a= [... 0] e     c= [... 0]
% N é o número de linhas
N = length(d);
 
% Modifica o primeiro coeficiente de cada linha
c(1) = c(1) / b(1); % Risco de divisão por zero.
d(1) = d(1) / b(1); 
 
for n = 2:1:N
    temp = b(n) - a(n) * c(n - 1);
    c(n) = c(n) / temp;
    d(n) = (d(n) - a(n) * d(n - 1)) / temp;
end
 
% Agora voltando a substituição.
x(N) = d(N);
for n = (N - 1):-1:1
    x(n) = d(n) - c(n) * x(n + 1);
end
end