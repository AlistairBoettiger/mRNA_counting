%%                       GB_simp_promoter.m
% 
% Alistair Boettiger                                   Date Begun: 02/03/11
% Levine Lab                                        Last Modified: 02/04/11

%% Description
% Simulate simple promoter bursting to fit to hb data on fraction of
% activated nuclei, total mRNA, activation time and degredation rates.  

% Nchs ... number of chemical species
% q ... number of chemical reactions
%    knu ... 1-by-q row vector; rate constants of the chemical reactions;
%            knu = k/nu^(order-1), where nu=volume and order=order of the chemical reaction
%    c ... Nchs-by-q; c(j,i) = coeficient at the j-th chemical species at the reactant side of the i-th reaction
%    p ... Nchs-by-q; p(j,i) = coeficient at the j-th chemical species at the product side of the i-th reaction
%    A0 ... Nchs-by-1 column vector;
%           A0(j) = number of molecules of j-th chemical species at the initial time t=0
%    tfin ... final time; if the simulation raches tfin it stops

%clear all;


% Species
% 1. promoter_closed
% 2. promoter_open
% 3. mRNA

% Reactions
% -> promoter_open
% -> promoter_closed
% promoter_closed -> promoter_open
% promoter_open -> promoter_closed
% promoter_open -> promoter_open + mRNA
% mRNA -> 0


% Mean burst size should be: b =  k_tx / k_off
% Inter-burst interval should be: t_off =  1/k_on
Nreal = 100;

tfin = 60*20;

% mu_M = (k_on*k_off/(k_off + k_on)*b )/k_d

mu_M = 1000; % This is MEASURED directly
t_off = 30; % lifetime of off state = time to switch on.  This is MEASURED  
f = .5; % fraction on.  THIS IS MEASURED
k_d = 1/(5*60);  % 5 min -- should be able to MEASURE this!

k_on = 1/t_off;  % defined by lifetime of off state
k_off = k_on*(1/f-1); % chosen so f fraction are on;

k_tx = k_d*mu_M*(k_off+k_on)/(k_on);
%k_tx = b*k_off;


b = k_tx/k_off;  disp(['b = ',num2str(b)]);

N = 100;
xx = linspace(0,1,N);

% p = gamma(k_on + k_off)/(gamma(k_on)*gamma(k_off))*(1/k_tx)*x.^(k_on-1).*(1-x).^(k_off-1);
% 
% figure(1); clf; plot(x,p);

k = [0 0 k_on k_off k_tx k_d];

c = [0, 0, 1, 0, 0, 0;
     0, 0, 0, 1, 1, 0;
     0, 0, 0, 0, 0, 1];
 
 
p = [1, 0, 0, 1,0,0;
     0, 1, 1, 0,1,0;
     0, 0, 0, 0,1,0];
A0 = [1;
      0;
      0];
  


%ode45

% plot Nreal realizations

figure(1); clf;

figure(2); clf;

dots = zeros(Nreal,1); 

for n = 1 : Nreal
  [A,t] = gillespiessa(k, c,p, A0,tfin);
  
%  [A,t] = GibsonBruck(k, c,p, A0,tfin);

 
%   figure(1);
%   plot( t,A(1,:), 'ro');
%   hold on;
%   plot( t,A(2,:), 'g.');
%    plot( t,A(3,:), 'b.');

   dots(n) = A(3,end);
end % for n

  figure(1);
  %plot( t,A(1,:), 'ro');
  hold on;
  plot( t,A(2,:), 'g.');
  plot( t,A(3,:), 'b.');

figure(2); clf; hist(dots);
cov = std(dots)/mean(dots);
title(['cov=',num2str(cov,4)]);
