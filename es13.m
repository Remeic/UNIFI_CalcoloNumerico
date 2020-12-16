% Esercizio: 13
% Scaricare la function cremat1 al sito:
% http://web.math.unifi.it/users/brugnano/appoggio/cremat1.m
% che crea sistemi lineari m × n, con m ? n, la cui soluzione (nel senso dei minimi
% quadrati) `e il vettore x = ( 1 ...  n)^(T)
% Eseguire, quindi, il seguente script Matlab
% per testare le function dei precedenti esercizi

format long 
i=1;
for n = 5:10
xx = [1:n]';
    for m = n:n+10
        [A,b] = cremat1(m,n);
        QR = myqr(A);
        x = qrsolve(QR,b);
        rst(i,1) = i;
        rst(i,4) = m;
        rst(i,2) = n;
        rst(i,3) = norm(x-xx);
        i=i+1;
    end
end

colNames = {'m','n','norma','norm_x_xx'};
tableResult = array2table(rst,...,
    'VariableNames',colNames);
disp(tableResult);


function [A,b] = cremat1(m,n)
%
%  [A,b] = cremat1(m,n)       Crea una matrice A nxn ed un termine noto b,
%                             in modo che la soluzione del sistema lineare,
%                             A*x=b, nel senso dei minimi quadrati, sia 
%                             x = [1,2,...,n]^T.
%
rng(0);
A = rand(m,n);
[q,r] = qr(A);
b = r*[1:n]';
b(n+1:m) = rand(m-n,1);
b = q*b;
return
end



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

