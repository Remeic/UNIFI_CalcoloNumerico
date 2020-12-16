% Esercizio 15
% Scrivere un programma che implementi efficientemente il calcolo del
% polinomio interpolante di Hermite su un insieme di ascisse distinte.


% Confrontare i codici degli esercizi 14-17 per approssimare la funzione
% f(x) = sin(x) sulle ascisse xi = i*pi/n, i = 0, 1, . . . , n, per n = 1, 2, . . . , 10. Graficare
% l'errore massimo di approssimazione verso n (in semilogy), calcolato su una griglia
% uniforme di 10001 punti nell'intervallo [0, pi].


risultatoHermite = []
x = linspace(0,pi,10001);
f = @(x) sin(x);
f1 = @(x) cos(x);
y = feval(f,x);


for n=2:10
   xi = (0:n) * pi/n;
    fi = feval(f, xi);
    ff1 = feval(f1, xi);
    hermiteY =  hermite(xi, fi, ff1, x)
    errorHermite = max(abs(y - hermiteY));
    risultatoHermite = [risultatoHermite errorHermite];
end
x=(2:10);
semilogy(x,risultatoHermite);



function y = hermite(xi, fi, fi1, x)
% y = hermite(xi, fi, fi1, x)
% Funzione che calcola il polinomio interpolante per una data coppia di dati
% (xi,fi) nei punti del vettore x, utilizzando il metodo di Hermite
% 
% Parametri:
%    xi: ascisse di interpolazione
%    fi: valori della funzione e la sua derivata nelle ascisse di interpolazione
%    x: ascisse dove valutare il polinomio, con il metodo di hermite
% Restituisce:
%    y: valutazione delle ascisse in x 
% 
   if (length(xi) - length(fi) ~= 0) | (length(xi) - length(fi1) ~= 0)
       error("La lunghezza dei vettori xi, fi e fi1 deve essere la medesima")
   end
   if isempty(x) 
      error("Il vettore x non deve essere vuoto") 
   end
   for u=1:length(x)-1
       if(x(u) == x(u+1)) 
           fprintf("Ascisse uguali in posizione %d - %d \n",x(u),x(u+1));
           error("Le ascisse devono essere distinte");
       end
   end
   xi = reshape([xi; xi], [], 1)';
   fi = reshape([fi; fi1], [], 1)';
   n = length(xi)-1;
   % Calcolo differenze divise
   for i=n:-2:3
       fi(i) = (fi(i)-fi(i-2))/((xi(i)-xi(i-1)));
   end
   for j=2:n
       for i=n+1:-1:j+1
           fi(i) = (fi(i)-fi(i-1))/((xi(i)-xi(i-j)));
       end
   end
   % Algoritmo di Horner
   y = fi(n+1)*ones(size(x));
   for k=1:length(x) 
       for i=n:-1:1
           y(k) = y(k)*(x(k)-xi(i))+fi(i);
       end
   end
   return 
end