
%%                          State-site model, full simulation
%
% Alistair Boettiger                                  Date Begun: 01/18/11
% Levine Lab                                        Last Modified: 01/18/11

function sim_all_binding
clear all;
fout = '/Users/alistair/Documents/Berkeley/Levine_Lab/Projects/Enhancer_Modeling/Results';


% States

N = 6; % number of binding sites


bound = cell(1,N); 
AllStates = [];
nstates = zeros(1,N);
for n = 1:N+1
    nstates(n) = nchoosek(N,n-1); 
    bound{n}  = zeros(1,nstates(n));
    AllStates = [AllStates, bound{n}];
end

AllStates = AllStates';
S = length(AllStates);

k0 = .03; % much better this is slow and later ones are fast.  
k1 = .1; 
k2= 1;
km1 = 100;
km2 = 25;
km3 = 1;




G = zeros(S);  % Transition Matrix

%        0    1      2     3     12      13     23     123
 Gu= [-k0*3   k0     k0    k0    0         0     0       0;  % 0
        0   -2*k1    0     0     k1       k1     0       0;  % 1    
        0    0   -2*k1     0     k1       0      k1      0;  % 2   
        0    0      0    -2*k1   0        k1     k1      0;  % 3   
        0    0      0      0    -k2       0      0      k2;  % 12   
        0    0      0      0     0       -k2     0      k2;  % 13   
        0    0      0      0     0        0     -k2     k2;  % 23   
        0    0      0      0     0        0      0      0];  % 123   
    

    
    %    0    1     2     3    12      13    23      123
 Gl= [   0    0     0    0      0       0     0      0;  % 0
        km1 -km1    0     0     0       0     0      0;  % 1    
        km1   0   -km1    0     0       0     0      0;  % 2   
        km1   0     0  -km1     0       0      0     0;  % 3   
        0    km2   km2    0  -2*km2   0       0      0;  % 12   
        0    km2    0    km2    0   -2*km2    0      0;  % 13   
        0     0    km2   km2    0       0    -2*km2  0;  % 23   
        0     0     0     0   km3     km3     km3     -3*km3];  % 123  

% Need an algorithm to build this correctly for arbitrary N.    
  
ks = (1:n+1)*.1;
km = fliplr( logspace(-2,2,n+1) ); 
[Gu, Gl] = nFillingMatrix(N,ks,km);

Cs = 60; % different concentration values to try  =  % number of nuclei in cc13 embryo  
AllStates(1) = 1;
state = zeros(Cs,S);
c = linspace(1,Cs,Cs);
%c = logspace(-4,2,Cs);


c_max = 58.5 ; % molecules per micron square at x = 0 
%  necessary for c = 4.8 at x = .5
a_eff = 0.003; % effective size of primary enhancer
tau = 60*3; % integration time seconds
D = 4.5; % um^2/s


pts = Cs; % number of nuclei in cc13 embryo
x = linspace(0,1,pts);
c = c_max*D*a_eff*tau*exp(-x*5);
figure(1); clf; plot(x,c,'g'); ylabel('molecules/um^3'); 

half_pt = round(pts/2); % the middle point, x = .5
Thresh = c(half_pt); %
p_see_above_thresh = zeros(1,pts);
p_see_below_thresh = zeros(1,pts);
for i = 1:pts
   p_see_above_thresh(i) = gammainc(c(i),Thresh,'lower');  % 1 - upper is lower  % I think this floor is artifactual.
    
   % p_see_above_thresh(i) = gammainc(c(i),floor(Thresh + .1*Thresh),'lower');  % 1 - upper is lower
    % p_see_below_thresh(i) = gammainc(c(i),floor(Thresh - .1*Thresh),'upper');  % 1 - upper is lower
end       


figure(2); clf; plot(x,p_see_above_thresh,'c.');
xlabel('cell position'); ylabel('probability of expression');
title(['a=',num2str(a_eff,3),'um tau = ',num2str(tau), 's  D=',num2str(D), 'c_{1/2} = 4.8 #/um^3']);
set(gcf,'color','w'); 

% hold on; plot(x,p_see_below_thresh,'r.');
% ect_on = sum(p_see_above_thresh(half_pt:end))/half_pt
% 
% 
% 
% ect_off = 1-sum(p_see_above_thresh(1:half_pt))/half_pt
% 
% ect_on = sum(p_see_above_thresh(half_pt+1:end))/half_pt
% 
% y = zeros(1,100);
% c = linspace(0,7,100);
% for i=1:100
%     y(i) = gammainc(c(i),3.8,'lower');
% end
% figure(2); clf; plot(c,y);

%% Theory:

Ke1 = km1/k0; Ke2 = km2/k1;  Ke3 = km3/k2;

p  = c.^3./(  (Ke1*Ke2*Ke3) + (Ke2*Ke3)*c + (Ke3)*c.^2 + c.^3);
p_opt =  c.^3./(  (Ke1*Ke2*Ke3) +  c.^3);
p_opt6 =  c.^6./(  (Ke1*Ke2*Ke3).^2 +  c.^6);

figure(1); hold on; plot(c,p,'k--','linewidth',3);
plot(c,p_opt,'r--','linewidth',3);
plot(c,p_opt6,'g--','linewidth',3);
  
%% Direct solution  
% for i=1:Cs
%     G =  (Gl+Gu*c(i));
%     state(i,:) = G'  \ (AllStates);  % solve for equilibrium occupancy
% end
% 

% %% Euler integration
% T = 1000;
% dt = .001;
% check = zeros(1,Cs); 
% checkU = zeros(1,Cs); 
% checkL = zeros(1,Cs); 
% mleak = zeros(1,Cs);
% 
% stateT = zeros(S,Cs,T);
% for i=1:Cs
%     G = (Gl+Gu*c(i));
%     check(i) = sum(sum(G,2));
%     checkU(i) = sum(sum(Gu*c(i),2)); % some round off error of order 10E-15
%     
%     temp = AllStates*c(i);
%     leak = zeros(1,T); 
%     for t=1:T
%          temp = temp + dt*G'*temp;  % solve for equilibrium occupancy
%         stateT(:,i,t) = temp;
%         leak(t) = sum(temp);
%     end
%     state(i,:) = stateT(:,i,T)';
%     mleak(i) = mean(leak);
% end
% figure(5); clf; colordef white; set(gcf,'color','w');
% plot(state);

%% ODE method
T = 100;

timedat = cell(2,Cs);
state = NaN*zeros(Cs,S);
for i=1:Cs
    G = (Gl+Gu*c(i))*.1;
    R0 = AllStates;
    
    [t R] = ode15s(@sys, [0 T], R0,[],G);
    save([fout,'/','test.mat']); 
    state(i,:) = R(end,:);
    timedat{1,i} = t;
    timedat{2,i} = R;
%     
%     figure(2); clf; plot(R); legend('0','1','2','3','12','13','23','123');
% hold on; plot(sum(R,2),'k--','linewidth',2)

%pause(.1);

end

%
figure(1); clf; colordef white; set(gcf,'color','w');

col = jet(S);

figure(1); clf;
for j=1:S
    plot(c,state(:,j),'.','color',col(j,:)); hold on;
end
plot(c,p,'k--');
plot(c,p_opt,'r--','linewidth',1);
plot(c,p_opt6,'g--','linewidth',1);
legend('0','1','2','3','12','13','23','123','det. theory','opt for N=3','opt for N=6');

figure(2); clf; plot(timedat{1,45},timedat{2,45});

figure(3); clf;
for j=1:S
    plot(x,state(:,j),'.','color',col(j,:)); hold on;
end
 plot(x,p,'k--');
plot(x,p_opt,'r--','linewidth',1);
plot(x,p_opt6,'g--','linewidth',1);
legend('0','1','2','3','12','13','23','123','det. theory','opt for N=3','opt for N=6');


           function dRdt = sys(t,R,G)  
                dRdt = G'*R;

% load([fout,/,'test.mat']); 

