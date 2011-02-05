% Chemical system with SNIPER bifurcation:
%       0 --> Y      with rate constant k1d
%       Y --> X      with rate constant k2d
%       X --> 0      with rate constant k3d
%      2X --> 3X     with rate constant k4d
%      3X --> 2X     with rate constant k5d
%   X + Y --> X + 2Y with rate constant k6d
%  2X + Y --> 2X     with rate constant k7d

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

fout = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Enhancer_Modeling/Data/';

k1d = 12;   % sec^(-1) mm^(-3)
k2d = 1;    % sec^(-1)
k3d = 33;   % sec^(-1)
k4d = 11;   % sec^(-1) mm^3
k5d = 1;    % sec^(-1) mm^6
k6d = 0.6;  % sec^(-1) mm^3
k7d = 0.13; % sec^(-1) mm^6
vol = 40; % mm^3
k1nu = k1d*vol;   % sec^(-1)
k2nu = k2d;       % sec^(-1)
k3nu = k3d;       % sec^(-1)
k4nu = k4d/vol;   % sec^(-1)
k5nu = k5d/vol^2; % sec^(-1)
k6nu = k6d/vol;   % sec^(-1)
k7nu = k7d/vol^2; % sec^(-1)

q = 2;
Nchs = 2;

k = [0 0 .1 .1 .1];

c = [0, 0, 1, 0, 0;
     0, 0, 0, 1, 1;
     0, 0, 0, 0, 0];
 
 
p = [1, 0, 0, 1,0;
     0, 1, 1, 0,1;
     0, 0, 0, 0,1];
A0 = [2;
      0;
      0];
  
tfin = 1000;

%ode45

% plot Nreal realizations
Nreal = 1;
figure(1)
clf;
%hold on;
figure(2)
clf;
%hold on;
for n = 1 : Nreal
  [A,t] = gillespiessa(k, c,p, A0,tfin);
  
  [A,t] = GibsonBruck(k, c,p, A0,tfin);

 
  figure(1);
  plot( t,A(1,:), 'ro');
  hold on;
  plot( t,A(2,:), 'g.');
   plot( t,A(3,:), 'b.');

  fprintf('done [%fsec].\n', toc);
end % for n


