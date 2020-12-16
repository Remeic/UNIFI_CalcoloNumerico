% Scrivere una function che implementi efficientemente il metodo delle
% potenze.

format longg
tol = 10^(-10);
for  i=1:1:50
    n=i*10;
    result(i,1) = n;
    [result(i,3) result(i,2)] = metodoPotenze(matrSparsa(n), tol, ones(n,1));
end


columnsName = {'n','Numero_Iterazioni','Stima_Autovalore'};
tableResults = array2table(result,...
    'VariableNames',columnsName);
disp(tableResults);

function [lambda, i] = metodoPotenze(A, tol, x0, maxit) 
% Utilizzo: [lambda, i] = metodoPotenze(A, [tol, [x0, [maxit]]])
% Restituisce l'autovalore dominante della matrice A e il numero di 
% iterazioni necessarie per calcolarlo.
%   
% Parametri:
%       A: matrice utilizzata per il calcolo
%       [tol]: tolleranza dell' approssimazione (specificata o omessa)
%       [x0]: vettore di partenza (specificato o omesso)
%       [maxit]: numero massimo di iterazioni
% Restituisce:
%       lambda: matrice quadrata nxn sparsa
%       i: numero di iterazioni
%
    [m,n] = size(A); % Inizializzo m ed n con la dimensione della matrice
    if m ~= n % Controllo se la matrice è una matrice quadrata
        error('Errore: La matrice deve essere quadrata.'); 
    end
    if nargin <= 1
        tol = 10^(-6);
    elseif (tol >= 0.1) || (tol <= 0)
            error("Errore: Tolleranza specificata non corretta");
    end
    if nargin <= 2 % Genero il vettore iniziale se non specificato
        x=rand(n,1);
    else
        x = x0;
    end
    if nargin <= 3 % Inizializzo il numero di iterazioni massime se non specificate
        maxit = ceil(-log(tol))*n;
    end
    x = x/norm(x); 
    lambda = inf; 
    for i = 1:maxit
        lambda0 = lambda;
        v = A*x; 
        lambda = x'*v;
        err = abs(lambda-lambda0); % Calcolo dell'errore 
        if err < tol*(1+abs(lambda)) % Se il valore dell'errore è accettabile interrompo il ciclo
            break;
        end
        x = v/norm(v);
    end
    if err > tol*(1+abs(lambda)) % Se il valore dell'errore non è accettabile dopo le iterazioni termino
        warning('la tolleranza richiesta non è stata raggiunta.');
    end
end

function A = matrSparsa(n)
% Utilizzo: A = matrSparsa(n)
% Genera la matrice quadrata sparsa nxn, con n maggiore di 10.
%   
% Parametri:
%       n: numero di righe/colonne della matrice quadrata sparsa
% Restituisce:
%       A: matrice quadrata nxn sparsa
%    
    if n < 10
        error('n deve essere maggiore di 10.');
    end
    
    d = ones(n,1)*4;
    A = spdiags(d,0,n,n);
    
    d = ones(n,1)*(-1);
    A = spdiags(d,1,A);
    A = spdiags(d,-1,A);
    
    if n >= 10
        A = spdiags(d,9,A);
        A = spdiags(d,-9,A);
    end
end