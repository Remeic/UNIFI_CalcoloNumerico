% Esercizio: 11
% Scrivere una function Matlab che, data in ingresso una matrice
% A appartenente R m×n con m >= n = rank(A), restituisca una matrice, QR, che contenga
% l'informazione sui fattori Q ed R della fattorizzazione QR di A:
% function QR = myqr(A)
% Curare particolarmente la scrittura e l?efficienza della function

function QR = myqr(A)
% Utilizzo: [A] = myqr(A)
% Calcola la fattorizzazione QR della matrice A.
% Parametri:
%    A: la matrice da fattorizzare
% Restituisce:
%    QR: una matrice QR scritta con la parte significativa di R e la parte
%       significativa dei vettori di Householder normalizzati con la prima   
%       componente unitaria
    [m,n]=size(A); % Inizializzo m ed n con la grandezza di A
    for i=1:n
    alpha = norm(A(i:m, i), 2); % Rank della matrice
    if alpha==0 
        error('La matrice non ha rank massimo');
    end
    if(A(i,i))>=0
        alpha = -alpha;
    end
    v1 = A(i,i)-alpha;
    A(i,i) = alpha;
    A(i+1:m, i) = A(i+1:m, i)/v1;
    beta = -v1/alpha;
    A(i:m, i+1:n) = A(i:m, i+1:n) -(beta*[1; A(i+1:m, i)])*([1 A(i+1:m, i)']*...
                    A(i:m,i+1:n));
    end
    QR = A;
end
