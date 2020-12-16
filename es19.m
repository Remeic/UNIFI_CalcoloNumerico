% Esercizio: 19
% Calcolare (numericamente) la costante di Lebesgue per i polinomi
% interpolanti di grado n = 2, 4, 6, . . . , 40, sia sulle ascisse equidistanti che su quelle di
% Chebyshev (utilizzare 10001 punti equispaziati per valutare la funzione di Lebesgue).
% Graficare convenientemente i risultati ottenuti. Spiegare, quindi, i risultati ottenuti
% approssimando la funzione
% f(x) = 1 + (1 + x^2) x fra -5 e 5
% utilizzando le ascisse equidistanti e di Chebyshev precedentemente menzionate 
% (tabulare il massimo errore valutato su una gliglia 10001 punti equidistanti nell?intervallo
% [?5, 5]).

format short e 


f=@(x) 1/(1+x^2);
z=linspace(-5, 5, 10001)';
lebesgueConstanti=zeros(20, 1);
%lebesgueEquidistanti(f,z,lebesgueConstanti);
lebesgueCheby(f,z,lebesgueConstanti);

function lebEqui = lebesgueEquidistanti(f,z,lebesgueConstanti)
% Utilizzo: lebesgueEquidistanti(f,z,lebesgueConstanti)
% Funzione che calcola la costante di lebesgue per polinomi
% interpolanti di grado n:2...40 su ascisse equidistanti
%
% Parametri: 
%   f: Funzione da approssimare
%   z: punti di valutazione
%   lebesgueConstanti: vettore per memorizzare le costanti calcolate
%
iteraz=0;
    for n=2:2:40
        xi = zeros (1 , n+1);
        fxi = zeros(1,n+1);
        for i = 1:n+1
            xi(i) = -5+(i-1)*(10/n);
        end
        for i = 1:n+1
            fxi(i) = feval(f, xi(i));
        end
        fx=zeros(1, 10001); %valori nelle 10001 ascisse uniformi
        lebesgueConstanti(n/2) = norm(lebesgue(xi),inf);
        for k=1:10001
            fx(k)=feval(f,z(k));
        end
        y=lagrange(xi, fxi,z);
        err=max(abs(fx-y'));
        iteraz = iteraz +1;
        rst(iteraz,1) = n;
        rst(iteraz,2)= err;
        rst(iteraz,3)=lebesgueConstanti(n/2);
    end
    x=2:2:40
    plot(x,lebesgueConstanti);
    title('Costante di Lebesgue (Ascisse Equidistanti)');
     colNames = {'n','errore','norma'};
tableResult = array2table(rst,...,
    'VariableNames',colNames);
disp(tableResult);
end

function lebCheb = lebesgueCheby(f,z,lebesgueConstanti)
% Utilizzo: lebesgueEquidistanti(f,z,lebesgueConstanti)
% Funzione che calcola la costante di lebesgue per polinomi
% interpolanti di grado n:2...40 su ascisse equidistanti
%
% Parametri: 
%   f: Funzione da approssimare
%   z: 
%   lebesgueConstanti: vettore per memorizzare le costanti calcolate
%
    iteraz = 0;
    for n=2:2:40
        xi = zeros (1 , n+1);
        fxi = zeros(1,n+1);
        for i = 1:n+1
            xi = chebyshev(n, -5, 5);
            xi = sort(xi);
        end
        for i = 1:n+1
            fxi(i) = feval(f, xi(i));
        end
        fx=zeros(1, 10001); %valori nelle 10001 ascisse uniformi
        lebesgueConstanti(n/2) = norm(lebesgue(xi),inf);
        for k=1:10001
            fx(k)=feval(f,z(k));
        end
        y=lagrange(xi, fxi,z);
        err=max(abs(fx-y'));
        iteraz = iteraz +1;
        rst(iteraz,1) = n;
        rst(iteraz,2)= err;
        rst(iteraz,3)=lebesgueConstanti(n/2);
    end
    x=2:2:40
    plot(x,lebesgueConstanti);
    title('Costante di Lebesgue (Ascisse di Chebyshev)');
    colNames = {'n','errore','norma'};
tableResult = array2table(rst,...,
    'VariableNames',colNames);
disp(tableResult);
end


function y = lagrange(xi, fi, x)
% y = lagrange(xi, fi, x)
% Funzione per il calcolo del polinomio interpolante per le coppie di dati (xi,fi) nei
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


function leb = lebesgue(puntiInterp)
% Utilizzo: leb = lebesgue(puntiInterp)
% Funzione per il calcolo della costante di lebesgue 
%
% Parametri:
% puntiInterp: punti di interpolazione
%
% Restituisce:
% y: Costante di lebesgue
       npunti = length(puntiInterp);
       lin=zeros(10001,1);    
       x=linspace(-5,5,10001);
       for j=1:10001
           k=0;
           for i=1:npunti
               val=abs(lagrangePol(puntiInterp,x(j),i));
               k=k+val;
           end
           lin(j,1)=k;
       end
           leb=lin;
return
end


function valPol = lagrangePol(z,x,i)
% Utilizzo: lagrangePol(z,x,i)
% Funzione per il calcolo dell' i-esimo polinomio di Lagrange
%
% Parametri:
%
%   z: punti di interpolazione
%   x: punto su cui effettuare la valutazione
%   i: indice del polinomio
%
% Restituisce:
%   valPol: Valore del polinomio in x

    n = length(z); m = length(x);
    valPol = prod(repmat(x,1,n-1)-repmat(z([1:i-1,i+1:n]),m,1),2)/...
    prod(z(i)-z([1:i-1,i+1:n])); 
return
end


function xi = chebyshev(n, a, b)
   % xi = ceby(n, a, b)
   % Calcola le ascisse di Chebyshev per il polinomio di grado n
   % trasformate in [a,b]
   %
   % Parametri:
   %    n: grado del polinomio interpolante
   %    a: estremo sinistro
   %    b: estremo destro
   %
   % Restituisce:
   %    xi: ascisse di Chebyshev
  
   if n<=0
       error('Il numero n inserito deve essere maggiore di 0.')
   end
   xi = cos((2*[0:n]+1)*pi/(2*n+2));
   xi = ((a+b)+(b-a)*xi)/2;
   return
end


