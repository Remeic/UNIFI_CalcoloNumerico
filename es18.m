
x=linspace (0 , pi , 10001);
fx=zeros (1 , 10001);
fx(1:10001)=sin(x(1:10001));
mv=zeros(10, 1);


for n=1:10
    xi = zeros(n+1,1);
    fi = zeros(n+1,1);
    fd = zeros(n+1,1);
    for i=0:10
        xi(i+1) = (i*pi)/n;
        fi(i+1) = sin(xi(i+1));
        fd(i+1) = cos(xi(i+1));
    end
    
    % Lagrange
%     y = lagrange(xi,fi,x);
%     maxErr = max(abs(fx-y));
%     mv(n) = maxErr;
    
    % Hermite
    y = sistemaSplineNaturale(xi,fi);
    maxErr = max(abs(fx-y));
    mv(n) = maxErr;
    display(y)
end

semilogy(mv);

function y = SplineNotaKnot( xi, fi )
% y = SplineNotaKnot( xi, fi, x )
% Calcola e valuta i valori della spline not-a-knot 
% relativa alle coppie di dati assegnati.
% Parametri:
%    xi: ascisse di interpolazione
%    fi: valori della funzione, valutati nelle ascisse xi
%    x: punti di valutazione della spline 
% Restituisce:
%    y: valutazione dei punti x calcolati sulla spline Not a Knot

        mi = sistemaSplineNotAKnot(xi, fi);        
        y = valutaSplineNotAKnot(xi, fi, mi, x);
    return
end

function result = differenzeDivise(x, f)
% Parametri:
% - f: funzione;
% - x: vettore delle ascisse.
% Restituisce:
% - result: l'ultima differenza divisa calcolata
%   essa include tutte le ascisse
    n = length(x);
    for i = 1:n-1
        for j = n:-1:i+1
            f(j)=(f(j)-f(j-1))/(x(j)-x(j-i));
        end
    end
    result = f(end);
end



function m = sistemaSplineNotAKnot(xi, fi)
% Calcola i coefficienti m da applicare all'espressione della spline cubica
% not-a-knot.
% Parametri:
%    xi: ascisse di interpolazione
%    fi: valori della funzione, valutati nelle ascisse xi
% Restituisce:
%    m: coefficienti necessari a calcolare l'espressione della spline
%    not-a-knot
    
    % calcolo delle phi, zi e diff_div
    % n:= grado del polinomio interpolante
    n = length(xi)-1;
    hi = xi(2) - xi(1);
    hi1 = xi(3) - xi(2);
    phi(1) = hi/(hi+hi1);
    zi(1) = hi1/(hi+hi1);
    diffDiv(1) = differenzeDivise(xi(1:3), fi(1:3));
    for i = 2:n-1 
        hi = xi(i) - xi(i-1);
        hi1 = xi(i+1) - xi(i);
        phi(i) = hi/(hi+hi1);
        zi(i) = hi1/(hi+hi1);
        diffDiv(i) = differenzeDivise(xi(i-1:i+1), fi(i-1:i+1));
    end
    diffDiv(i+1) = differenzeDivise(xi(i:i+2), fi(i:i+2));
    % espando il vettore diffDiv per avere lunghezza n+1
    diffDiv = [6*diffDiv(1), 6*diffDiv, 6*diffDiv(end)];
    % i = 1
    u(1) = 1;
    w(1) = 0;
    l(1) = phi(1) / u(1);
    u(2) = 2 - phi(1);
    w(2) = zi(1) - phi(1);
    l(2) = phi(2) / u(2);
    u(3) = 2 - l(2)*w(2);
    % i = 4:n-1
    for i = 4:n-1
        w(i-1) = zi(i-2);
        l(i-1) = phi(i-1) / u(i-1);
        u(i) = 2 - l(i-1)*w(i-1);
    end
    % i = n
    w(n-1) = zi(n-2);
    l(n) = (phi(n-1) ...
        - zi(n-1)) / ...
        u(n-1);
    u(n) = 2 - zi(n-1) - l(n)*w(n-1);
    % i = n+1
    w(n) = zi(n-1);
    l(n+1) = 0;
    u(n+1) = 1;
    % 1) Ly = 6*diffDiv
    y(1) = diffDiv(1);
    for i = 2:n+1
        y(i) = diffDiv(i) - l(i-1)*y(i-1);
    end
    % 2) Um = y
    m = zeros(n+1, 1); 
    m(n+1) = y(n+1) / u(n+1);
    for i = n:-1:1
        m(i) = (y(i) - w(i) * m(i+1)) / u(i);
    end
    m(1) = m(1) - m(2) - m(3);
    m(end)= m(n+1)-m(n)-m(n-1);
    return
end

function y = valutaSplineNotAKnot(xi, fi, mi, x)
% Valuta i punti della spline cubica not-a-knot nei punti assegnati x
% Parametri:
%    xi: ascisse di interpolazione
%    fi: valori della funzione, valutati nelle ascisse xi
%    mi: coefficienti necessari a calcolare l'espressione della spline
%    x: punti in cui valutare la spline
% Restituisce:
%    y: valutazioni nei punti x della spline.
    
    if xi(1)>x(1) || xi(end)<x(end)
        error('I punti di valutazione non rientrano nel dominio della spline.');
    end
    if length(xi)~=length(fi)
        error('I vettori xi e fi devono avere la stessa lunghezza.');
    end
    % raccolgo tutti i sottoinsiemi di punti da valutare con le relative
    % functions, i punti x vegono ordinati
    sort(x);
    y = zeros(length(x),1);    
    lastIndex = 1;
    k = 2;
    for j = 1 : length(x)
        if x(j) <= xi(k-1)
            j = j + 1;
        else
            if j ~= lastIndex
                % calcolo della spline relativa al k-1° intervallo
                hi = xi(k) - xi(k-1);
                ri = fi(k-1) - ((hi^2)/6) * mi(k-1);
                qi = (fi(k) - fi(k-1))/hi - (hi/6) * ...
                    (mi(k) - mi(k-1));
                spline = @(x) (((x - xi(k-1)).^3) .* mi(k) + ((xi(k) - x).^3) .* mi(k-1))./(6*hi)...
                  + qi .* (x - xi(k-1)) + ri;
                % calcolo le valutazioni della spline
                y(lastIndex : j-1) = spline(x(lastIndex : j-1));
                lastIndex = j;
            end
            k = k+1;
        end
    end
    if j ~= lastIndex
        hi = xi(end) - xi(end-1);
        ri = fi(end-1) - ((hi^2)/6) * mi(end-1);
        qi = (fi(end) - fi(end-1))/hi - (hi/6) * (-mi(end-1));
        spline = @(x) (((x - xi(end-1)).^3) .* mi(end) + ((xi(end) - x).^3) .* mi(end-1))./(6*hi)...
                  + qi .* (x - xi(end-1)) + ri;
        y(lastIndex:j-1) = spline(x(lastIndex:j-1));
    end
end




function m = sistemaSplineNaturale( xi , fi )
% INPUT:
%    xi: ascisse di interpolazione
%    fi: valori della funzione, valutati nelle ascisse xi
% OUTPUT:
%    m: coefficienti necessari a calcolare l'espressione della spline
%    cubica

    % calcolo delle phi, zi e diff_div (lunghi n-1)
    % n:= grado del polinomio interpolante
    n = length(xi)-1;
    phi = zeros(n-2, 1);
    zi = zeros(n-2, 1);
    d = zeros(n-1, 1);
    hi = xi(2) - xi(1);
    hi1 = xi(3) - xi(2);
    phi(1) = hi/(hi+hi1);
    zi(1) = hi1/(hi+hi1);
    d(1) = differenzeDivise(xi(1:3), fi(1:3));
    for i = 2:n-2 
        hi = xi(i) - xi(i-1);
        hi1 = xi(i+1) - xi(i);
        phi(i) = hi/(hi+hi1);
        zi(i) = hi1/(hi+hi1);
        d(i) = differenzeDivise(xi(i-1:i+1), fi(i-1:i+1));
    end
    d(i+1) = differenzeDivise(xi(i:i+2), fi(i:i+2));
    % calcolo delle li ed ui che completano il sistema lineare della spline
    % naturale
    u(1) = 2;
    for i = 2:n-1 
        l(i) = phi(i-1)/u(i-1);
        u(i) = 2-l(i)*zi(i-1);
    end
    % risoluzione del sistema lineare LU per trovare gli m
    % 1) L*y = 6*b
    % il vettore d delle diffDiv viene sovrascritto
    d(1) = 6*d(1);
    d(2:n-1) = 6*d(2:n-1)-l(2:n-1)*d(1:n-2);
    % 2) U*m = y
    m(n-1) = d(n-1)/u(n-1);
    m(n-2:-1:1) = (d(n-2:-1:1)-zi(n-2:-1:1) * m(n-1:-1:2))/u(n-2:-1:1);
end

function y = valutaSplineNaturale(xi, fi, mi, x)
    % INPUT:
    %    xi: ascisse di interpolazione
    %    fi: valori della funzione, valutati nelle ascisse xi
    %    mi: coefficienti necessari a calcolare l'espressione della spline
    %    x: punti in cui valutare la spline
    % OUTPUT:
    %    y: valutazioni nei punti x della spline.
    
    if xi(1)>x(1) || xi(end)<x(end)
        error('I punti di valutazione non rientrano nel dominio della spline.');
    end
    if length(xi)~=length(fi)
        error('I vettori xi e fi devono avere la stessa lunghezza.');
    end
    % raccolgo tutti i sottoinsiemi di punti da valutare con le relative
    % functions, i punti x vegono ordinati
    mi = [0  mi  0];
    sort(x);
    y = zeros(length(x),1);    
    lastIndex = 1;
    k = 2;
    for j = 1 : length(x)
        if x(j) <= xi(k-1)
            j = j + 1;
        else
            if j ~= lastIndex
                % calcolo la spline relativa al k-1° intervallo
                hi = xi(k) - xi(k-1);
                ri = fi(k-1) - ((hi^2)/6) * mi(k-1);
                qi = (fi(k) - fi(k-1))/hi - (hi/6) * ...
                    (mi(k) - mi(k-1));
                spline = @(x) (((x - xi(k-1)).^3) .* mi(k) + ((xi(k) - x).^3) .* mi(k-1))./(6*hi)...
                  + qi .* (x - xi(k-1)) + ri;
                % calcolo le valutazioni della spline
                y(lastIndex : j-1) = spline(x(lastIndex : j-1));
                lastIndex = j;
            end
            k = k+1;
        end
    end
    % valutazione degli ultimi punti
    if j ~= lastIndex
        hi = xi(end) - xi(end-1);
        ri = fi(end-1) - ((hi^2)/6) * mi(end-1);
        qi = (fi(end) - fi(end-1))/hi - (hi/6) * (-mi(end-1));
        spline = @(x) (((x - xi(end-1)).^3) .* mi(end) + ((xi(end) - x).^3) .* mi(end-1))./(6*hi)...
                  + qi .* (x - xi(end-1)) + ri;
        y(lastIndex:j-1) = spline(x(lastIndex:j-1));
    end
    
end


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
   xi = reshape([xi; xi], [], 1)'; % tipo trasforma il vettore in matrice ????
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