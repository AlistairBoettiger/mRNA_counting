

% Here is the implementation of Gibson-Bruck algorithm.
% 
% I based my implementation upon the implementation of the Gillespie algorithm of Tomas Vejchodsky
% see: www.math.cas.cz/vejchod/gillespiessa.html
% 
% Some of my own comments are in French, post back for eventual explanations.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,t,prop] = GibsonBruck(knu,c,p,A0,tfin)

%

% Next Reaction Method Stochastic Simulation Algorithm for a system of q
% chemical reactions

% implémenté selon le modèle gillespiessa

% Nchs ... number of chemical species
% q ... number of chemical reactions
%
% i-th reaction with the rate constant knu(i), i=1,2,...,q:
% c(1,i) A(1) + c(2,i) A(2) + ... + c(Nchs,i) A(Nchs)
% --->
% p(1,i) A(1) + p(2,i) A(2) + ... + p(Nchs,i) A(Nchs)
%
% INPUT:
% knu ... 1-by-q row vector; rate constants of the chemical reactions;
% knu = k/nu^(order-1), where nu=volume and order=order of the chemical reaction
% c ... Nchs-by-q; c(j,i) = coeficient at the j-th chemical species at the reactant side of the i-th reaction
% p ... Nchs-by-q; p(j,i) = coeficient at the j-th chemical species at the product side of the i-th reaction
% A0 ... Nchs-by-1 column vector;
% A0(j) = number of molecules of j-th chemical species at the initial time t=0
% tfin ... final time; if the simulation raches tfin it stops
%
% OUTPUT:
% A ... Nchs-by-N; A(j,m) = number of j-molecules in time interval [t(m), t(m+1)]
% A(:,1) == A0(:);
% t ... 1-by-N;
% prop (ajouté par moi) propensity pour chaque état


% NRM SSA
% =============
% cf feuille "Algorithme de Gillespie appliqué au modéle de LV stochastique"

% Initial checks
if (size(A0,2) > 1)
  error('The initial numbers of molecules A0 must form a COLUMN vector');
end
Nchs = size(A0,1);

if (size(knu,1) > 1)
  error('The rate constants knu must form a ROW vector');
end
q = size(knu,2);

if (size(c) ~= [Nchs,q])
  error('The reactant coeficients must be Nchs-by-q but Nchs=%d, q=%d and size(c)=(%d,%d)', Nchs,q,size(c,1),size(c,2));
end

if (size(p) ~= [Nchs,q])
  error('The product coeficients must be Nchs-by-q but Nchs=%d, q=%d and size(p)=(%d,%d)', Nchs,q,size(p,1),size(p,2));
end

% assume maximum of Nalloc steps of the algorithm
% allocate memory for Nalloc steps
Nalloc = 10000;

A = zeros(Nchs,Nalloc);
t = zeros(1,Nalloc);
prop = zeros(q,Nalloc);% chaque ligne correpond à une réaction, on y écrit les propensions associées
tau = zeros(q,Nalloc); % putative times (un par réaction)

% queue: vecteur avec q lignes (# reactions chimiques) et 2 colonnes
% (indice de la réaction et moment de la réaction)
queue = zeros(q,2);

% set the initial values
m = 1;
A(:,1) = A0;
t(1) = 0;

% propension (# reactions)
alpha = zeros(1,q);

% vecteur aléatoire uniforme (# reactions)
r = rand(1,q);
%pause

while t(m) <= tfin
    
    if m > 1
        
        %x=x
        % espèces touchées par la réaction x
        esp_x = find( (c(:,x)>0)|(p(:,x)>0) );%
        % quelles sont les reactions où ces espèces apparaissent
        react_x = find(sum(c(esp_x,:)>0,1)>0);%
        
        
        %alpha = zeros(1,q);
        % propensions actualisées (propension qui chg => dépend de x)
        for i = react_x % react touchées
            %i = i
            alpha(i) = knu(i);%
            for j=1:Nchs % go thru all chemical species
                %j=j%
                for s=1:c(j,i) % ... alphai=knu*A*(A-1)*(A-2)*...*(A-c(j,i)+1)
                    alpha(i) = alpha(i) * (A(j,m)-s+1);%
                end % for s
            end % for j
            prop(i,m) = alpha(i);%
            %pause
        end % for i
        
% % affichage
% m = m
% propen = alpha
% queue = queue
% temps = t(1:m)
%
% pause
        
        
        
        react_x = react_x(logical(react_x ~=x));%
        % temps de réactions actualisés
        for i = react_x % react touchées
            %i = i
            %alpha_old_ = alpha_old(i)
            %alpha_ = alpha(i)
            %tau_i = tau(i)
            %t_m = t(m-1)
            queue(i,2) = alpha_old(i)/alpha(i)*(queue(i,2)-t(m)) + t(m);%
            %pause
        end % for i
        
        % temps de réactions pour x
        rr = rand(1);
        queue(x,2) = -log(rr)/alpha(x) + t(m);%
        
        % sauvegarde propension
        alpha_old = alpha;
        
% % affichage
% m=m
% propen = alpha
% queue = queue
% temps = t(1:m)
%
% pause
        
       
        
    else % tp t=1, tirer au hasard pour toutes les propensions
        
        % propensions, tirées de gillespiessa
        alpha = zeros(1,q);
        for i=1:q % go thru all reactions and compute the propensity function alpha(i)
            alpha(i) = knu(i);
            for j=1:Nchs % go thru all chemical species
                for s=1:c(j,i) % ... alphai=knu*A*(A-1)*(A-2)*...*(A-c(j,i)+1)
                    alpha(i) = alpha(i) * (A(j,m)-s+1);
                end % for s
            end % for j
            prop(i,m) = alpha(i);
        end % for i
        
        % sauvegarde propension
        alpha_old = alpha;
        
        % temps de réactions
        tau(:,m) = -(log(r)./alpha)';%
        
        % queue
        queue = [(1:q)',tau(:,m)];%
        
% % affichage
% m = m
% r = r
% alpha = alpha
% tau_1 = tau(:,1)
% queue = queue
% temps = t(1:m)
% pause
    end
    
    % prochain état
    queue_triee = sortrows(queue,2);%
    x = queue_triee(1,1);% % l'indice de la reaction qui a lieu est celui pr
                          % lequel le temps de séjour est minimal
    A(:,m+1) = A(:,m) - c(:,x) + p(:,x);%
    t(m+1) = queue_triee(1,2);
    prop(:,m+1) = prop(:,m);%
    
% if mod(m+1,1000)==0
% mm = m+1
% tautau = tau(x)
% tt = t(m+1)
% AA = A(:,m+1)
% pause
% end
    
% % affichage
% m=m
% queue_triee = queue_triee
% A_m = A(:,m)
% A_mp1 = A(:,m+1)
% t_m = t(m)
% t_mp1 = t(m+1)
% prop_m = prop(:,m)
% prop_mp1 = prop(:,m+1)
% temps = t(1:(m+1))
% pause
 
        
    % incrément
    m = m + 1;%
    
    
    % memory management
    if (m >= Nalloc)
        %fprintf('Doubling memory for data. (Reallocation.) t(%d)=%f', m, t(m));
        %tic
        Aaux = zeros(Nchs,2*Nalloc);
        taux = zeros(1,2*Nalloc);
        Aaux(:,1:Nalloc) = A;
        taux(1:Nalloc) = t;
        A = Aaux;
        t = taux;
        clear Aaux taux;
        Nalloc = 2*Nalloc;
        %fprintf(' done. [%fsec]\n', toc);
    end
end

% cutting the zeros at the end of arrays
A = A(:,1:m);
t = t(1:m);
prop = prop(:,1:(m-1));

% postprocessing
if (t(m) > tfin)
  t(m) = tfin;
end
for j=1:Nchs
   if (A(j,m) < 0)
        A(j,m) = 0;
   end
end 