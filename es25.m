% Esercizio: 25
clear
clc
format longg
tol = 10^(-10);
for  i=1:1:50
    n=i*10;
    [result(i,1) result(i,2)] = metodoPotenze(matrSparsa(n), tol, ones(n,1));
end
semilogy(result(:,2));
%semilogy(result(:,1));
columnsName = {'Stima_Autovalore_Dominante','Numero_Iterazioni'};
tableResults = array2table(result,...
    'VariableNames',columnsName);
disp(tableResults);


function [lambda, i] = metodoPotenze(A, tol, x0, maxit) 
% Utilizzo: [lambda, i] = metodoPotenze(A, [tol, [x0, [maxit]]])
% Funzione che restituisce l'autovalore di modulo massimo della matrice A e il numero di 
% iterazioni necessarie per il suo calcolo.
% Parametri:
%   A: matrice utilizzata per il calcolo
%   [tol]: tolleranza dell' approssimazione (specificata o omessa)
%   [x0]: vettore di partenza (specificato o omesso)
%   [maxit]: numero massimo di iterazioni
% Restituisce:
%    lambda: matrice quadrata nxn sparsa
%    i: numero di iterazioni eseguite
%
    [m,n] = size(A); %Inizializzo m ed n con la dimensione della matrice
    if m ~= n  % Controllo se la matrice è una matrice quadrata
        error('La matrice deve essere quadrata.'); 
    end
    if nargin <= 1 % Inizializzo la tolleranza se non è specificata
        tol = 10^(-6);
    elseif (tol >= 0.1) || (tol <= 0) % Controllo di robustezza sulla tolleranza
        error("Tolleranza specificata non corretta");
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
    lambda = inf; % Inizializzo lambda con il valore infinito
    for i = 1:maxit
        lambda0 = lambda;
        v = A*x; 
        lambda = x'*v;
        err = abs(lambda-lambda0); % Calcolo dell'errore 
        if err < tol*(1+abs(lambda))% Se il valore dell'errore è accettabile interrompo il ciclo
            break;
        end
        x = v/norm(v);
    end
    if err > tol*(1+abs(lambda))  % Valore dell'errore non accettabile dopo le iterazioni
        warning('la tolleranza richiesta non è stata raggiunta.');
    end
end

function A = matrSparsa(n)
% A = matrSparsa(n)
% Genera la matrice quadrata sparsa nxn, con n maggiore di 10.
%   
% Parametri:
%    n: numero di righe/colonne della matrice quadrata sparsa
% Restituisce:
%    A: matrice quadrata nxn sparsa
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