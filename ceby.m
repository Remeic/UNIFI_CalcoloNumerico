function xi = ceby(n, a, b)
    % xi = ceby(n, a, b)
    % Calcola le ascisse di Chebyshev per il polinomio di grado n
    % trasformate in [a,b]
    xi = cos((2*[0:n]+1)*pi/(2*n+2));
    xi = ((a+b)+(b-a)*xi)/2;
    return
end

function xi = ceby(n, a, b)
    % xi = ceby(n, a, b)
    % Calcola le ascisse di Chebyshev per il polinomio di grado n
    % trasformate in [a,b]
    xi = cos((2*[0:n]+1)*pi/(2*n+2));
    xi = ((a+b)+(b-a)*xi)/2;
    return
end

function y = lagrange(xi, fi, x)
   % y = lagrange(xi, fi, x)
   % Calcolo del polinomio interpolante per le coppie di dati (xi,fi) nei
   % punti del vettore x
    
   % Controllo lunghezza di xi e fi
   if length(xi) - length(fi) ~= 0
       return ;
   end
   
   n = length(xi)-1;
   y = zeros(size(x));
   % per ogni ascissa in x
   for i=1:length(x)
       % calcolo del valore del polinomio
       for k=1:n+1
           lkn = 1;
           % calcolo di lkn
           for j=1:n+1
               if (j~=k)
                   lkn = lkn * ((x(i)-xi(j))/(xi(k)-xi(j)));
               end
           end
           y(i) = y(i) + fi(k)*lkn;
       end
   end
       
   return 
end