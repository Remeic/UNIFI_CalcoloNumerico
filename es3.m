% Esercizio: 3
% Nel primo caso abbiamo prima la moltiplicazione seguita dalla divisione, in questo 
% contesto non si presentano errori poichè non u non assume un valore minore del 
% valore minimo rappresentabile da matlab ( numero strano ), differente è il secondo caso
% in cui effettuando prima le divisioni si arriva ad una condizione di underflow 
% e il numero viene approssimato.
% Matlab in questo caso usa un sistema di recovery ( gradual underflow ) che diminuisce i bit
% della mantissa a favore dell'esponente
% L'errore generato dal primo ciclo di divisioni viene propagato dal
% successivo ciclo di moltiplicazioni

format long e
n=75;
u=1e-300;

rst(1,1)=u;
for i=1:n
    u=u*2;
end
rst(1,2)=u;
for i=1:n
    u=u/2;
end
rst(1,3)=u;
colNames = {'Iterazione','Moltiplicazione','Divisione'};
tableResult = array2table(rst,...,
    'VariableNames',colNames);
disp(tableResult);


rst(1,1)=u;
for i=1:n
    u=u/2;
end
rst(1,2)=u;
for i=1:n
    u=u*2;
end
rst(1,3)=u;
colNames = {'Iterazione','Moltiplicazione','Divisione'};
tableResult = array2table(rst,...,
    'VariableNames',colNames);
disp(tableResult);


