% Esercizio: 23
% Scrivere una function Matlab che, data in ingresso la matrice QR
% creata dalla function del precedente esercizio, ed il termine noto del sistema lineare
% Ax = b, ne calcoli la soluzione nel senso dei minimi quadrati:
% function x = qrsolve(QR,b)
% Curare particolarmente la scrittura e l?efficienza della function

function x = qrsolve(QR,b)
    % x = qrsolve(QR,b)
    % Calcola la soluzione di Ax=b con A matrice QR
    % Parametri:
    %   A: Matrice QR
    %   b: vettore colonna
    % Restituisce:
    %  b: soluzione del sistema lineare sovradeterminato
    [m,n] = size(QR);
    qtrasp=eye(m);
    for i=1:n
        qtrasp = [eye(i-1) zeros(i-1, m-i+1); zeros(i-1, m-i+1)' ...  
                 (eye(m-i+1)-(2/norm([1; QR(i+1:m, i)], 2)^2)*([1; QR(i+1:m, i)]...
                 *[1 QR(i+1:m, i)']))]*qtrasp;
    end
    x = triangSup(triu(QR(1:n, :)), qtrasp(1:n,:)*b);
end



function [b] = triangSup(A, b)
    % Utilizzo: [b] = diagonale(A, b)
    % Calcola la soluzione di Ax=b con A matrice triangolare superiore
    % Input: 
    %   A: Matrice triangolare superiore
    %   b: vettore colonna
    % Output:
    %  b: soluzione del sistema
    for i=length(A):-1:1
        for j=i+1:length(A)
            b(i)=b(i)-A(i,j)*b(j);
        end
        if A(i,i)==0
            error('Matrice singolare')
        else
            b(i)=b(i)/A(i,i);
        end
    end
end