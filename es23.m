% Esercizio: 23
% Sapendo che If=0atan(30)(1 + tan2 (x))dx
% tabulare il numero dei punti richiesti dalle formule adattative dei trapezi 
% e di Simpson per approssimare If, utilizzate con tolleranze 
% tol= 10-i i=2,...,8 assieme ai relativi errori.

f = @(x) (1+(tan(x).^2));

If = trapad(0,atan(30),f,10^(-8))
If2 = simpad(0,atan(30),f,10^(-8))


function If = trapad( a, b, f, tol, fa, fb )
% Utilizzo: If = trapad( a, b, f, tol)
% Calcola ricorsivamente l'integrale della funzione, nell'intervallo prescelto, 
% usando la formula dei trapezi adattiva.
%
% Input: 
%   a: estremo sinistro
%   b: estremo destro
%   f: funzione integranda
%   tol: tolleranza prefissata
% Output:
%  If: l’approssimazione dell’integrale definito della funzione

    % Controlli di robustezza:
    % - a deve essere minore di b
    if a>=b
        error("Intervallo di integrazione non corretto.")
    end
    global count;
    if nargin<=4
        fa = f(a);
        fb = f(b);
        count = 2;
    end
    h = b-a;
    x1 = (a+b)/2;
    f1 = f(x1);
    count = count + 1;
    I1 = (h/2)*(fa+fb);
    If = (I1+h*f1)/2;
    e = abs(If-I1)/3;
    if e>tol
       If =  trapad( a, x1, f, tol/2, fa, f1) + trapad( x1, b, f, tol/2, f1, fb);
    end
    assignin('base',"nPoint",count)
end



function If = simpad( a, b, f, tol, fa, f1, fb )
% If = simpad( a, b, f, tol)
% Funzione che calcola ricorsivamente l'integrale della funzione data, 
% in un intervallo specificato, utilizzando la formula adattiva di Simpson
%
% Parametri:
%   a: estremo sinistro dell?intervallo
%   b: estremo destro dell?intervallo
%   f: funzione integranda
%   tol: tolleranza prefissata ( Obbligatoriamente diversa da 0 )
% Restituisce:
%  If: Approssimazione dell?integrale definito della funzione in un dato 
%  intervallo
   
   if tol <= 0, error("Tolleranza specificata non corretta"), end
   if a >= b, error("Intervallo di integrazione non corretto."), end
   x1 = (a + b) / 2; % Calcolo del punto medio x1
   global count;
   if nargin <= 4 % Prima iterazione non ricorsiva
       fa = f(a); % Valutazione della funzione nei punti a,b,x1
       fb = f(b);
       f1 = f(x1);
	 count = 2; % Inizializzo count a 2 per i punti passati come arg.
    end
    h = (b - a) / 6;
    x2 = (a + x1) / 2;
    x3 = (x1 + b) / 2;
    count = count + 3; % Aggiungo a count il calcolo dei punti x1,x2,x3
    f2 = f(x2);
    f3 = f(x3);
    I1 = h*(fa+4*f1+fb); % Formula di Simpson
    If = .5*h*(fa+4*f2+2*f1+4*f3+fb); % Formula di Simpson Adattiva
    e = abs(If-I1)/15; % Calcolo dell?errore
    if e>tol % Chiamata ricorsiva alla funzione Simpad con nuovi intervalli
       If1 = simpad( a, x1, f, tol/2, fa, f2, f1);
       If2 = simpad( x1, b, f, tol/2, f1, f3, fb);
       If = If1 + If2;
    end
    assignin('base',"nPoint",count);  
end


