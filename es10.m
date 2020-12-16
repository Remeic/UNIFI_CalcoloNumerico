% Esercizio: 10
% Scaricare la function cremat al sito:
% http://web.math.unifi.it/users/brugnano/appoggio/cremat.m
% che crea sistemi lineari n × n la cui soluzione è il vettore x = ( 1 ...
% n) ^T

format short e
n = 10;
x = zeros(n,5);
for i = 1:5 
    [A,b] = cremat(n,i+10);
    [LU,p] = palu(A);
    x(:,i) = lusolve(LU,p,b);
end

colNames = {'i_1','i_2','i_3','i_4','i_5'};
tableResult = array2table(x,...,
    'VariableNames',colNames);
disp(tableResult);


function [A,b] = cremat(n,k,simme)
%
%  [A,b] = cremat(n,k,simme)  Crea una matrice A nxn ed un termine noto b,
%                             in modo che la soluzione del sistema lineare 
%                             A*x=b sia x = [1,2,...,n]^T.
%                             k non ve lo dico a cosa serve.
%                             simme, se specificato, crea una matrice
%                             simmetrica e definita positiva.
%
    if nargin==1
        sigma = 1/n; 
    else    
        sigma = 10^(-k);
    end    
    rng(0);
    [q1,r1] = qr(rand(n));
    if nargin==3
        q2 = q1';
    else    
        [q2,r1] = qr(rand(n));
    end    
    A = q1*diag([sigma 2/n:1/n:1])*q2;
    x = [1:n]';
    b = A*x;
    return
end



function [LU, p] = palu(A)
% Utilizzo: [LU, p] = palu(A)
% Funzione che calcola la fattorizzazione LU della matrice A con
% pivoting parziale di A

% Parametri:
%   A: la matrice da fattorizzare.
% Restituisce:
%   LU: la matrice fattorizzata LU;
%   p: vettore di permutazione

    n=size(A,1);
    p=(1:n);
    for i=1:n-1
        [mi, ki] = max(abs(A(i:n, i)));
        if mi==0 % Controllo il caso in cui la matrice sia singolare
            error('Matrice singolare');
        end
        ki = ki+i-1;
        if ki>i
            A([i ki], :) = A([ki i], :);
            p([i ki]) = p([ki i]);
        end
        A(i+1:n, i) = A(i+1:n, i)/A(i, i);
        A(i+1:n, i+1:n) = A(i+1:n, i+1:n) -A(i+1:n, i)*A(i, i+1:n);
    end
    LU = A;
end


function b = lusolve(LU,p,b)
% x = lusolve(A,b,p)
% Funzione che calcola la soluzione di Ax=b con A matrice LU con la tecnica del pivoting
% Parametri:
%   A: Matrice matrice LU con pivoting (generata dalla funzione palu)
%   b: vettore colonna
% Restituisce:
%   b: soluzione del sistema
    P = zeros(length(LU));
    for i = 1:length(LU)
        P(i, p(i)) = 1;
    end
    b = triangInfUnit(tril(LU,-1)+eye(length(LU)), P*b);
    b = triangSup(triu(LU), b);
   
end


function b = triangInfUnit(A, b)
% Utilizzo: b = tri_inf_unit(A, b)
% Calcola la soluzione del sistema Ax=b con A matrice triangolare inferiore a
% diagonale unitaria.
% Parametri: 
%   A: Matrice triangolare inferiore
%   b: vettore dei termini noti
% Restituisce:
%  b: soluzione del sistema
   
    [n,m] = size(A);
    if n~=m
        error('Matrice non quadrata.')
    end
    if ~istril(A)
        error('Matrice non triangolare inferiore')
    end
    if ~all(diag(A)==1)
        error('Matrice non a diagonale unitaria')
    end
    for i=1:n
        for j=1:i-1
            b(i) = b(i) - A(i,j)*b(j);
        end
        if A(i,i)==0  % Controllo se la matrice è singolare
            error('Matrice singolare')
        else
            b(i) = b(i)/A(i,i);
        end
    end
end

function [b] = triangSup(A, b)
% Utilizzo: [b] = diagonale(A, b)
% Calcola la soluzione di Ax=b con A matrice triangolare superiore
% Parametri: 
%   A: Matrice triangolare superiore
%   b: vettore colonna
% Restituisce:
%  b: soluzione del sistema
    for i=length(A):-1:1
        for j=i+1:length(A)
            b(i)=b(i)-A(i,j)*b(j);
        end
        if A(i,i)==0 % Controllo se la matrice è singolare
            error('Matrice singolare')
        else
            b(i)=b(i)/A(i,i);
        end
    end
end

