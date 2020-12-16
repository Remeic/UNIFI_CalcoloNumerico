% Esercizio: 26-27
format longg
i=1;
tol = 10^(-10);
for n=10:10:500
    rst(i,1) = n;
    [x,k] = splittingGenerico(matrSparsa(n), ones(n,1), @MsolveGauss, tol, zeros(n,1));
    rst(i,2) = k;
    i = i+ 1;
end
plot(rst(:,1),rst(:,2),'-');
colNames = {'n','numero_iterazioni',};
tableResult = array2table(rst,...
    'VariableNames',colNames);
disp(tableResult); 


function [x,i,nr] = splittingGenerico(A, b, Msolve, tol, x0, maxit)
% [x,i] = jacobi(A, b, [tol, [xo, [maxit]]])
% Restituisce la soluzione del sistema lineare Ax=b approssimata con il 
% metodo utilizzato dalla funzione Msolve e il numero di iterazioni eseguite.
% Parametri:
%    A: matrice utilizzata per il calcolo
%    b: vettore dei termini noti
%    Msolve: funzione che implementa il metodo di risoluzione
%    [tol]: tolleranza dell' approssimazione (specificata o omessa)
%    [x0]: vettore di partenza (specificato o omesso)
%    [maxit]: numero massimo di iterazioni (specificato o omesso)
% Restituisce:
%    x: soluzione approssimata del sistema
%
    D = diag(A);
    if ~all(D) % Controllo la diagona principale non sia nulla
        error('La diagonale di A non deve avere elementi nulli');
    end
    n = length(b); %Inizializzo n con la lunghezza del vettore dei termini noti
    if nargin <= 3
        tol = 10^(-6);
    elseif (tol >= 0.1) || (tol <= 0)
            error("Tolleranza specificata non corretta");
    end
    if nargin <= 4 % Genero il vettore iniziale se non specificato
        x=rand(n,1);
    else
        x = x0;
    end
    if nargin <= 5 % Inizializzo il numero di iterazioni massime se non specificate
        maxit = ceil(-log(tol))*n;
    end
    
    for i = 1:maxit
        r = A*x - b;
        err = norm(r,inf);
        nr(i) = err;
        if err<=tol % Se l'errore è inferiore alla tolleranza termino il ciclo
            break;
        end
        r = Msolve(A,r); % Risolvo il sistema con il metodo implementato dalla Msolve
        x = x-r;
        
    end
    if err>tol % Se l
        warning('La tolleranza richiesta non è stata raggiunta.');
    end
end


function A = matrSparsa(n)
% A = matrSparsa(n)
% Genera la matrice quadrata sparsa nxn, con n maggiore eo uguale di 10.
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
        A = spdiags(d,10,A);
        A = spdiags(d,-10,A);
    end
end



function y = MsolveGauss(M, r)
    % y = MsolveGauss(M, r)
    % Restituisce la soluzione del sistema lineare Mx=r.
    % Ogni iterazione risolve un sistema triangolare inferiore
    y=r;
    n = length(r);
    for i = 1:n
       y(i) = y(i)/M(i,i);
       y(i+1 : n) = y(i+1 : n) - M(i+1 : n,i)*y(i);
    end
end


function y = MsolveJacobi(M, r)
    % y = MsolveJacobi(M, r)
    % Restituisce la soluzione del sistema lineare Mx=r.
    % Ogni iterazione risolve un sistema diagonale
    y = r./diag(M);
end
