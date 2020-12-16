% Esercizio: 22
% Scrivere due functions che implementino efficientemente le formule
% adattattive dei trapezi e di Simpson.

function If = trapad( a, b, f, tol, fa, fb )
% If = trapad( a, b, f, tol)
% Funzione che calcola ricorsivamente l'integrale della funzione data, 
% in un intervallo specificato, utilizzando la formula dei Trapezi adattiva
%
% Parametri: 
%   a: estremo sinistro dell'intervallo
%   b: estremo destro dell'intervallo
%   f: funzione integranda
%   tol: tolleranza prefissata ( Obbligatoriamente diversa da 0 )
% Restituisce:
%   If: Approssimazione dell?integrale definito della funzione in un dato 
%   intervallo
      if tol >= 0, error("Tolleranza specificata non corretta"), end
      if a >= b, error("Intervallo di integrazione non corretto."), end
      if nargin<=4 % Prima iterazione non ricorsiva
        fa = f(a); % Valutazione della funziona in a e b
        fb = f(b);
      end
      h = b-a; 
      x1 = (a+b)/2; % Calcolo del punto medio fra a e b
      f1 = f(x1);
      I1 = (h/2)*(fa+fb); % Formula dei Trapezi
      If = (I1+h*f1)/2; % Formula dei Trapezi adattiva
      e = abs(If-I1)/3;
      if e>tol
        If =  trapad( a, x1, f, tol/2, fa, f1) ...
        + trapad( x1, b, f, tol/2, f1,fb);
      end    
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
   
   if tol >= 0, error("Tolleranza specificata non corretta"), end
   if a >= b, error("Intervallo di integrazione non corretto."), end
   x1 = (a + b) / 2; % Calcolo del punto medio x1
   if nargin <= 4 % Prima iterazione non ricorsiva
       fa = f(a); % Valutazione della funzione nei punti a,b,x1
       fb = f(b);
       f1 = f(x1);
    end
    h = (b - a) / 6;
    x2 = (a + x1) / 2;
    x3 = (x1 + b) / 2;
    f2 = f(x2);
    f3 = f(x3);
    I1 = h*(fa+4*f1+fb); % Formula di Simpson
    If = .5*h*(fa+4*f2+2*f1+4*f3+fb); % Formula di Simpson Adattiva
    e = abs(If-I1)/15; % Calcolo dell?errore
    if e>tol % Chiamata ricorsiva alla funzione Simpad con nuovi intervalli
       If = simpad( a, x1, f, tol/2, fa, f2, f1) ...
           + simpad( x1, b, f, tol/2, f1, f3, fb);
    end
end