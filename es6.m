% Esercizio: 6
% Utilizzare le function del precedente esercizio per determinare una
% approssimazione della radice della funzione
% f(x) = x-e^(-x)cos(x/100)
% per tol = 10^(-i) i = 1, 2, . . . , 12, partendo da x0 = ?1. 
% Per il metodo di bisezione, utilizzare [-1,1], come intervallo di confidenza iniziale. 
% Tabulare i risultati, in modo da confrontare le iterazioni richieste da ciascun metodo. 
% Commentare il relativo costo computazionale.

f = @(x) x-(exp(-x))*cos(x/100);
% Derivata prima
derivate_f = @(x) (1/100)*exp(-x)*sin(x/100)+exp(-x)*cos(1/100)+1;

x0 = -1;
itmax = 30000;
rst = [0 0 0 0 0 ];
% Richiamo newton, corde e secanti per valori decrescenti di tolx
for i = 1:12
    tol = 10^-i;
    rst(i, 1) = tol;
    a = -1;
    b = 1;
    iterazioni = metodoBisezione(f,a,b,tol);
    rst(i,2) = iterazioni;
    iterazioni = metodoNewton(f, derivate_f, x0, tol, itmax);
    rst(i, 3) = iterazioni;
    [iterazioni,x] = metodoCorde(f, derivate_f, x0, tol, itmax);
    rst(i, 4) = iterazioni;
    iterazioni = metodoSecanti(f, derivate_f, x0, tol, itmax);
    rst(i, 5) = iterazioni;
end

colNames = {'tolx','Bisezione','Newton','Corde','Secanti'};
tableResult = array2table(rst,...,
    'VariableNames',colNames);
disp(tableResult);

function [iterazioni,x] = metodoBisezione(f, a, b, tolx)
% Utilizzo: [iterazioni,x] = metodoBisezione(f, a, b, tolx, itmax)
% Calcola dell'approssimazione di una radice della funzione
% METODO: Bisezione
%
% Parametri:
%   f: Funzione utilizzata
%   a: Estremo sinistro dell'intervallo di condifenza iniziale
%   b: Estremo destro dell'intervallo di condifenza iniziale
%   tolx: Tolleranza prefissata
% Restituisce:
%   iterazioni: numero di iterazioni eseguite 
%   x: Radice approssimata
    x = inf;
    fa = feval(f,a);
    fb = feval(f,b);
    if fa*fb > 0
        error("L' intervallo non contiene alcuna radice");
    end
    imax = ceil(log2(b-a)-log2(tolx));
    for i = 2:imax
        xi = (a+b)/2;
        fx = feval(f,x);
        if abs(xi-x) <= tolx*(1+abs(xi))
            break;
        elseif fa*fx < 0
            b = xi;
            fb = fx;
        else
            a = xi;
            fa = fx;
        end
        x = xi;
    end
    iterazioni = i;
end

function [iterazioni,x] = metodoNewton(f, f_der, x0, tolx, itmax)
% Utilizzo: [iterazioni,x] = metodoNewton(f, f_der, x0, tolx, itmax)
% Calcola dell'approssimazione di una radice della funzione
% METODO: Newton
%
% Parametri:
%   f: Funzione utilizzata
%   f_der: Derivata della funzione designata
%   x0: Punto di innesco
%   tolx: Tolleranza prefissata
%   itmax: Numero di iterazioni massime 
% Restituisce:
%   iterazioni: numero di iterazioni eseguite ( -1 nel caso di una non
%   convergenza)
%   x: Radice approssimata

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
    if abs(x-x0) <= tolx * (1 + x0) % Convergenza ottenuta
        iterazioni = i;
    else  % Convergenza non ottenuta
        warning("Il metodo (Newton) non converge");
        iterazioni = -1;
    end
end

function [iterazioni,x] = metodoCorde(f, f_der, x0, tolx, itmax)
% Utilizzo: [iterazioni,x] = metodoCorde(f, f_der, x0, tolx, itmax)
% Calcola dell'approssimazione di una radice della funzione
% METODO: Corde
%
% Parametri:
%   f: Funzione utilizzata
%   f_der: Derivata della funzione designata
%   x0: Punto di innesco
%   tolx: Tolleranza prefissata
%   itmax: Numero di iterazioni massime 
% Restituisce:
%   iterazioni: numero di iterazioni eseguite ( -1 nel caso di una non
%   convergenza)
%   x: Radice approssimata

    fx = feval(f, x0);
    f_derx = feval(f_der, x0);
    x = x0 - fx/f_derx;
    i = 0;
    while (i<itmax) && (abs(x-x0)>tolx*(1+x0))
        i = i+1;
        x0 = x;
        fx = feval(f, x0);
        x = x0 - fx/f_derx;
    end
    iterazioni = i;
    if abs(x-x0)> tolx*(1+x0) % Convergenza non ottenuta
        warning("Il metodo (Corde) non converge");
        iterazioni = -1;
    end
end

function iterazioni = metodoSecanti(f, f_der, x0, tolx, itmax)
    fx = feval(f, x0);
    f_derx = feval(f_der, x0);
    x = x0 - fx/f_derx;
    i = 0;
    while (i<itmax) && (abs(x-x0)>tolx*(1+x0))
        i = i+1;
        fx0 = fx;
        fx = feval(f, x);
        x1 = (fx*x0-fx0*x)/(fx-fx0);
        x0 = x;
        x = x1;
    end
    
    if abs(x-x0) <= tolx*(1+x0)  % Convergenza ottenuta
        iterazioni = i;
    else % Convergenza non ottenuta
        warning("Il metodo (Secanti) non converge");
        iterazioni = -1;
    end
end