% Esercizio: 21
% Uno strumento di misura ha una accuratezza di 10?6
% (in opportune
% unità di misura). I dati misurati nelle posizioni xi sono dati da yi come descritto
% dalla seguente tabella: ( vedi i vettori sotto )
% Calcolare il grado minimo, ed i
% relativi coefficienti, del polinomio che meglio approssima i precedenti dati nel senso
% dei minimi quadrati con una adeguata accuratezza. Graficare convenientemente i
% risultati ottenuti.

clear
clc
xi=[0.010;0.098;0.127;0.278;0.547;0.632;0.815;0.906;0.913;0.958;0.965];%vettore delle posizioni xi
yi=[1.003626;1.025686;1.029512;1.029130;0.994781;0.990156;1.016687;1.057382;1.061462;1.091263;1.096476];%misurazioni
n=length(xi);
tol=10^(-6);
g=calcolaGradoMinimo(n,xi,yi,tol);
V=vandermondeMatrix(n,g,xi);
QR=myqr(V);
a=qrsolve(QR,yi);
x1=linspace(0,1);
plot(xi,yi);

function gradoMinimo = calcolaGradoMinimo(n,xi,yi,tol)
% Utilizzo: gradoMinimo = calcolaGradoMinimo(n,xi,yi,tol)
% Calcola il grado minimo del polinomio nel senso dei minimi quadrati
%
% Parametri:
%   n: Grado della matrice Vandermonde n*n
%   xi: Ascisse per la creazione della matrice di Vandermonde
%   yi: Valore assunto dalla 
%   tol: tolleranza prefissata
%
% Restituisce: 
%   gradoMinimo: grado minimo calcolato nel senso dei minimi quadrati

    for gradoMinimo=1:10
        V=vandermondeMatrix(n,gradoMinimo,xi);
        QR=myqr(V);
        a=qrsolve(QR,yi);
        if(norm(V*a-yi,2)<tol)
            break;
        end
    end
    gradoMinimo=gradoMinimo-1;
    return
end


function VanderMatrix = vandermondeMatrix(n,m,xi)
% Utilizzo: VanderMatrix = vandermonde(m,xi) 
% Funzione per la creazione della matrice di vandermonde V 
% partendo da un vettore xi 
%
% Parametri:
%   n: numero righe di V
%   m: numero colonne di V
%   xi: elementi noti per la costruzione  della matrice
%
% Resistuisce:
%   V: matrice di Vandermonde
    for i=1:n-1
        for j=i+1:n
            if(xi(i)==xi(j))
                error("Errore: Le ascisse non sono tutte distinte");
            end
        end
    end
    for i=1:n
        for j=1:m
            z=j-1;
            VanderMatrix(i,j)=xi(i)^z;
        end
    end
    return
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
