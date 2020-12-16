% Esercizio: 8
% Scrivere una function Matlab che, data in ingresso una matrice A,
% restituisca una matrice, LU, che contenga l?informazione sui suoi fattori L ed U,
% ed un vettore p contenente la relativa permutazione, della fattorizzazione LU con
% pivoting parziale di A:

% function [LU,p] = palu(A)


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
            error('Errore: Matrice singolare');
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


function [b]= sistemaLUpivot(A,b,p)
    % [b]= sistemaLUpivot(A,b,p)
    % Calcola la soluzione di Ax=b con A matrice LU con pivoting
    % Parametri:
    % 	 A: Matrice matrice LU con pivoting (generata da palu)
    %   b: vettore colonna
    % Restituisce:
    %   b: soluzione del sistema
    P = zeros(length(A));
    for i = 1:length(A)
        P(i, p(i)) = 1;
    end
    b = tri_inf_unit(tril(A,-1)+eye(length(A)), P*b);
    b = triangSup(triu(A), b);
end
