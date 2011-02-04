clear all;


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

% Nchns different chemical species
mRNA = 1;
promoter_open = 0;
promoter_closed = 0; 


% Q reactions
% promoter_closed -> promoter_open
% promoter_open -> promoter_closed
% promoter_open -> mRNA + promoter
% promoter_closed + mRNA -> promoter_closed + mRNA  % dummy rxn


A0 = [mRNA;
      promoter_open;
      promoter_closed];

k = [.01,.01,1,1];

% reactions  
c = [1, 1, 0, 0;  % 1 promoter_closed
     0, 0, 0, 1;   
     0, 0, 1, 0];
 
 
 
 % products
 p = [0,0,1,1;
      1,0,0,0;
      0,1,0,0];
  

  T = 100; 
  
  [A,t,prop] = GibsonBruck(k,c,p,A0,T);
  
  figure(1); clf; 
  plot(t,A(1,:),'ko','linewidth',3);
  hold on; 
      plot(t,A(2,:),'gx','linewidth',3);  
      plot(t,A(3,:),'r.','linewidth',3);