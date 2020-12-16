% Esercizio: 7
% Calcolare la molteplicit`a della radice nulla della funzione
% f(x) = x^(2)sin(x^2)
% Confrontare, quindi, i metodi di Newon, Newton modificato, e di Aitken, 
% per approssimarla per gli stessi valori di tol del precedente esercizio (ed utilizzando il medesimo
% criterio di arresto), partendo da x0 = 1. Tabulare e commentare i risultati ottenut


clc
f = @(x) x.^2 * sin(x.^2);
f_der = @(x) 2*x.*(sin(x.^2) + x.^2*cos(x.^2));
x0=1;
itmax = 200;
molteplicita = 4;
for i=1:12   
    rst(i,1) = 10^(-i);
    [it,x] = newton(f,f_der,x0,10^(-i),itmax);
    rst(i,2) = it;
    rst(i,3) = x
    [it,x] = newton_mod(f,f_der,x0,molteplicita,10^(-i),itmax);
    rst(i,4) = it;
    rst(i,5) = x
    [it,x] = aitken(f,f_der,x0,10^(-i),itmax)
    rst(i,6) = it;
    rst(i,7) = x
end

colNames = {'Tol','It_Newton','x_Newton','It_NewtonMod','x_NewtonMod','It_Aitken','x_Aitken'};
tableResult = array2table(rst,...,
    'VariableNames',colNames);
disp(tableResult);



function [iterazioni,x] = newton(f, f_der, x0, tolx, itmax)
% Utilizzo [iterazioni,x] = newton_mod(f, f_der, x0, tolx, itmax)
% Calcola dell'approssimazione di una radice della funzione
% Metodo: Newton
% Parametri:
%   f: funzione utilizzata
%   f_der: derivata della funzione utilizzata f
%   x0: punto di innesco
%   molteplicità: molteplicità nota della radice
%   tolx: tolleranza prefissata
%   itmax: numero massimo di iterazioni
% Restituisce:
%   i: numero di iterazioni eseguite (-1 nel caso non converga)
%   x: radice approssimata    
    fx = feval(f, x0);
    f_derx = feval(f_der, x0);
    x = x0 - fx/f_derx;
    i = 0;
    while (i<itmax) && (abs(x-x0)>tolx)
        i = i+1;
        x0 = x;
        fx = feval(f, x0);
        f_derx = feval(f_der, x0);
        x = x0 - fx/f_derx;
    end
    if abs(x-x0)<=tolx
        iterazioni = i;
    else
        % Non abbiamo raggiunto la convergenza entro itmax iterazioni
        iterazioni = -1;
    end
end


function [iterazioni,x] = newton_mod(f, f_der, x0, molteplicita, tolx, itmax)
% Utilizzo [iterazioni,x] = newton_mod(f, f_der, x0, tolx, itmax)
% Calcola dell'approssimazione di una radice della funzione
% Metodo: Newton Modificato
% Parametri:
%   f: funzione utilizzata
%   f_der: derivata della funzione utilizzata f
%   x0: punto di innesco
%   molteplicità: molteplicità nota della radice
%   tolx: tolleranza prefissata
%   itmax: numero massimo di iterazioni
% Restituisce:
%   i: numero di iterazioni eseguite (-1 nel caso non converga)
%   x: radice approssimata    
    fx = feval(f, x0);
    f_derx = feval(f_der, x0);
    x = x0 - fx/f_derx;
    i = 0;
    while (i<itmax) && (abs(x-x0)>tolx)
        i = i+1;
        x0 = x;
        fx = feval(f, x0);
        f_derx = feval(f_der, x0);
        x = x0 - molteplicita*fx/f_derx;
    end
    if abs(x-x0)>tolx
        % Non abbiamo raggiunto la convergenza entro itmax iterazioni
        i = -1;
    end
    iterazioni = i;
end

function [iterazioni,x] = aitken(f, f_der, x0, tolx, itmax)
% Utilizzo [iterazioni,x] = aitken(f, f_der, x0, tolx, itmax)
% Calcola dell'approssimazione di una radice della funzione
% Metodo: Aitken
% Parametri:
%   f: funzione utilizzata
%   f_der: derivata della funzione utilizzata f
%   x0: punto di innesco
%   molteplicità: molteplicità nota della radice
%   tolx: tolleranza prefissata
%   itmax: numero massimo di iterazioni
% Restituisce:
%   i: numero di iterazioni eseguite (-1 nel caso non converga)
%   x: radice approssimata 
    i = 0;
    x = x0;
    vai = 1;
    while (i<itmax) && vai
        i = i + 1;
        x0 = x;
        fx = feval(f, x0);
        f_derx = feval(f_der, x0);
        x1 = x0 - fx/f_derx;
        fx = feval(f, x1);
        f_derx = feval(f_der, x1);
        x = x1 - fx/f_derx;
        x = (x*x0-x1^2)/(x-2*x1+x0);
        vai = abs(x-x0)>tolx;
    end
    if ~vai
        iterazioni = i;
    else
        % Non abbiamo raggiunto la convergenza entro itmax iterazioni
        iterazioni = -1;
    end
end