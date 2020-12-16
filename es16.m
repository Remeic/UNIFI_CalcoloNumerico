% Esercizio 16
% Scrivere un programma che implementi efficientemente il calcolo di
% una spline cubica naturale interpolante su una partizione assegnata.

% Confrontare i codici degli esercizi 14-17 per approssimare la funzione
% f(x) = sin(x) sulle ascisse xi = i*pi/n, i = 0, 1, . . . , n, per n = 1, 2, . . . , 10. Graficare
% l'errore massimo di approssimazione verso n (in semilogy), calcolato su una griglia
% uniforme di 10001 punti nell'intervallo [0, pi].

clear
clc
format  long e
risultatoSpline = [];
a = 0;
b = pi;
x = linspace(a,b,10001);
f = @(x) sin(x);
y = feval(f,x);

for n=2:10
    xi = (0:n) * pi/n;
    fi = feval(f, xi);
    splineY = splineCubicaNaturale(xi, fi, x);
    errorSpline = max(abs(y - splineY));
    risultatoSpline = [risultatoSpline errorSpline];
end
x=(2:10);
semilogy(x,risultatoSpline);

function y = splineCubicaNaturale(xi, fi, x)
% Utilizzo: splineNaturale = naturalSpline(xi, fi, x)
% Funzione che calcola e valuta i valori della spline cubica naturale 
% relativa alle coppie di dati assegnati.
%
% Parametri:
%   xi: punti di interpolazione
%   fi: valori della funzione, valutati nelle ascisse xi
%   x: punti di valutazione della spline
%
% Restituisce:
%   splineNaturale: valori di interpolazione della spline cubica naturale
%
    
    n = length(xi);
    splineNaturale = zeros(n, 1)';
    widthI = xi(2 : n) - xi(1 : n - 1);
    subDiag = (widthI(1 : end - 1)) ./ (widthI(1 : end - 1) + widthI(2 : end));
    superDiag = (widthI(2 : end)) ./ (widthI(1 : end - 1) + widthI(2 : end));
    divdiff = (fi(2 : n) - fi(1 : n - 1)) ./ widthI;
    divdiff = 6 * ((divdiff(2 : end) - divdiff(1 : end - 1)) ./ (xi(3 : end) - xi(1 : end - 2)));
    m = sistemaSplineNaturale(subDiag, superDiag, divdiff);
    [primaConst, secondaConst] = costantiIntegrazione(m, xi,fi, widthI);
    k=2;
    for j = 1 : length(x)
        for i = 2 : length(xi)
            if x(j) <= xi(i)
                k = i;
                break;
            end
        end
        y(j) = (((x(j) - xi(k - 1)) ^ 3) * m(k) + ...
          ((xi(k) - x(j)) ^ 3) * m(k - 1)) / ...
          (6 * widthI(k - 1)) + secondaConst(k - 1) * ...
          (x(j) - xi(k - 1)) + primaConst(k - 1);
    end
    return
end


function m = sistemaSplineNaturale(subDiag, superDiag, divDiff)
    % Utilizzo: m = sistemaSplineNaturale(subDiagCoeff, superDiagCoeff, divdiff)
    % Funzione che ritorna i coefficienti necessari per il calcolo della spline
    % cubica
    %
    % Parametri:
    %   subDiag: coeffiencenti della sotto diagonale della matrice 
    %   superDiag: coeffiencenti della sopra diagonale della matri
    %   diffDiv: differenze divise
    %
    % Restituisce:
    %   m: coefficienti necessari a calcolare l'espressione della spline
    %    cubica
    %
    n = length(superDiag) + 1;
    u(1) = 2;
    l = zeros(1, n - 2);
    m = zeros(1, n - 1);
    for i = 2 : n - 1
        l(i) = subDiag(i) / u(i - 1);
        u(i) = 2 - l(i) * superDiag(i - 1);
    end
    f(1) = divDiff(1);
    for i = 2:n - 1
        f(i) = divDiff(i) - l(i) * f(i - 1);
    end
    m(n - 1) = f(n - 1) / u(n - 1);
    for j = n - 2 : - 1 : 1
        m(j) = (f(j) - superDiag(j + 1) * m(j + 1)) / u(j);
    end
    m = [0 m 0];
    return
end

function [primaConst, secondaConst] = costantiIntegrazione(m, xi, fi, valI)
    % Utilizzo[ri, qi] = costantiIntegrazione(m, fi, xi, widthI)
    % Funzione che ritorna le costanti di integrazione della spline.
    %
    % Parametri:
    %   m: fattore m calcolato
    %   fi: valori della funzione, valutati nelle ascisse xi
    %   xi: ascisse di interpolazione
    %   valI: valore dell'i-esimo intervallo
    %
    % Restituisce:
    %   primaConst: valore della costante prima integrazione
    %   secondaConst: valore della costane seconda integrazione

    n = length(xi);
    primaConst = zeros(1, n-1);
    secondaConst = primaConst;
    for i = 2 : n
        primaConst(i - 1) = fi(i - 1) - (valI(i - 1) ^ 2) / 6 * m(i - 1);
        secondaConst(i - 1) = (fi(i) - fi(i - 1)) / ...
                    valI(i - 1) - valI(i - 1) / 6 * (m(i) - m(i - 1));
    end
    return
end



