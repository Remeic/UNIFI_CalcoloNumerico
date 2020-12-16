% Esercizio: 9
% Scrivere una function Matlab che, data in ingresso la matrice LU ed il
% vettore p creati dalla function del precedente esercizio (palu), ed il termine noto del sistema
% lineare Ax = b, ne calcoli la soluzione:
% function x = lusolve(LU,p,b)

function b = lusolve(LU,p,b)
    % x = lusolve(A,b,p)
    % Calcola la soluzione di Ax=b con A matrice LU con la tecnica del pivoting
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

