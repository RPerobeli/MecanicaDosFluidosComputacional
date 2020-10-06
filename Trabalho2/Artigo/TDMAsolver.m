function x = TDMAsolver(a,b,c,d)
% a, b, c s�o os vetores colunas para a compress�o da matriz tridiagonal, d � o vetor � direita
% a= [... 0] e     c= [... 0]
% N � o n�mero de linhas
N = length(d);
 
% Modifica o primeiro coeficiente de cada linha
c(1) = c(1) / b(1); % Risco de divis�o por zero.
d(1) = d(1) / b(1); 
 
for n = 2:1:N
    temp = b(n) - a(n) * c(n - 1);
    c(n) = c(n) / temp;
    d(n) = (d(n) - a(n) * d(n - 1)) / temp;
end
 
% Agora voltando a substitui��o.
x(N) = d(N);
for n = (N - 1):-1:1
    x(n) = d(n) - c(n) * x(n + 1);
end
end