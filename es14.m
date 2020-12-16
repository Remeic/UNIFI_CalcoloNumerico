% Esercizio 14
% Scrivere un programma che implementi efficientemente il calcolo del
% polinomio interpolante su un insieme di ascisse distinte.

% Confrontare i codici degli esercizi 14-17 per approssimare la funzione
% f(x) = sin(x) sulle ascisse xi = i*pi/n, i = 0, 1, . . . , n, per n = 1, 2, . . . , 10. Graficare
% l'errore massimo di approssimazione verso n (in semilogy), calcolato su una griglia
% uniforme di 10001 punti nell'intervallo [0, pi].

risultatoLagrange = []
x = linspace(0,pi,10001);
f = @(x) sin(x);
y = feval(f,x);


for n=2:10
   xi = (0:n) * pi/n;
    fi = feval(f, xi);
    lagrangeY =  lagrange(xi, fi, x)
    errorHermite = max(abs(y - lagrangeY));
    risultatoLagrange = [risultatoLagrange errorHermite];
end
x=(2:10);
semilogy(x,risultatoLagrange);


function y = lagrange(xi, fi, x)
% y = lagrange(xi, fi, x)
% Calcolo del polinomio interpolante per le coppie di dati (xi,fi) nei
% punti del vettore x, con il metodo di lagrange.
% Metodo utilizzato: Lagrange.
%
% Parametri:
% xi: ascisse di interpolazione
% fi: valori della funzione nelle ascisse di interpolazione
% x: ascisse dove valutare il polinomio, Metodo: lagrange
% Restituisce:
% y: valutazione delle ascisse in x 
%
   if length(xi) - length(fi) ~= 0
       error("La lunghezza dei vettori xi e fi non é la stessa.")
   end
   if isempty(x) 
      error("Vettore x vuoto.") 
   end
   n = length(xi)-1;
   for i=1:n % Controllo che tutte le ascisse siano distinte
        for j=i+1:n
            if xi(i)==xi(j) 
                error("Le ascisse non sono tutte distinte");
            end
       end
   end
   y = zeros(size(x));
   % per ogni ascissa in x
   for i=1:length(x)
       % calcolo del valore del polinomio
       for k=1:n+1
           lkn = 1;
           % calcolo di lkn
           for j=1:n+1
               if (j~=k)
                   lkn = lkn.*((x(i)-xi(j))/(xi(k)-xi(j)));
               end
           end
           y(i) = y(i) + fi(k)*lkn;
       end
   end
   return 
end