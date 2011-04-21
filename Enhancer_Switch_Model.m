%%                      Simplest Enhancer Switching Model

% Alistair Boettiger                                   Date Begun: 03/19/11
% Levine Lab                                        Last Modified: 04/20/11

%%

%% Actually interesting: (?)

% 10% off, 10% off, 1% off.
% imagine 1 chance to turn on, on cells make 200 +/- sqrt(200) transcripts
% off celsl make 0.

FAC = .9; % fraction of active cells 
mpc = 20; % mean mRNA per cell
Ncells = 250; % number of cells; 

sigma = std([mpc*ones(1,round(Ncells*FAC)) + sqrt(mpc)*(.5-rand(1,round(Ncells*FAC)))  , 0*ones(1,round(Ncells*(1-FAC) ))]);
mu = mean([mpc*ones(1,round(Ncells*FAC)) + sqrt(mpc)*(.5-rand(1,round(Ncells*FAC))), 0*ones(1,round(Ncells*(1-FAC)))]);
cov_1x = sigma/mu 

FAC = 1-(1-FAC)^2;
sigma = std([mpc*ones(1,round(Ncells*FAC)) + sqrt(mpc)*(.5-rand(1,round(Ncells*FAC))) , 0*ones(1,round(Ncells*(1-FAC)))]);
mu = mean([mpc*ones(1,round(Ncells*FAC)) + sqrt(mpc)*(.5-rand(1,round(Ncells*FAC))) , 0*ones(1,round(Ncells*(1-FAC)))]);
cov_2x = sigma/mu 

cov_1x/cov_2x

%%
%mean_dots =  120; %225;%   300;
% 5; % 2; % average time open (in seconds)
t_loop = 20;
t_close = 30; 
reload_rate = 1/4; % 16; % average number of polII reloaded per second
dots = zeros(1,500);  
time = 0; i = 0;
aveT_open =  t_close/t_loop; 
t_life = 6*60; 

Tot_time = 20*60;
T_op = Tot_time*t_close/(t_close +t_loop);
 mean_dots = reload_rate*T_op*exp(-(Tot_time)/t_life); %*(1-t_loop/Tot_time) %  t_close*
%(1-t_loop/Tot_time)


% Tot_time = mean_dots/reload_rate;

mean_burst = aveT_open*reload_rate;
mean_openings = Tot_time/aveT_open;

for n=1:500
    time = 0; 
    while time < Tot_time;
        T_loop = -log(rand)*t_loop;
        time = time + T_loop;
        T_open = -log(rand)*t_close;
        dots(n) = dots(n) + sum(cumsum(-log(rand(1,1000))/reload_rate )<T_open);
        dots(n) = dots(n)*exp(-(T_loop+T_open)/t_life); 
        time = time + T_open;
    end
end
figure(3); clf; hist(dots);
Mr = mean(dots);
Sr = std(dots);
CoV = std(dots)/mean(dots);


disp(['Mean Burst Size = ',  num2str(mean_burst,3) '  Mean promoter openings = ',num2str(mean_openings,3)  ]); 

disp(['Mean = ',num2str(Mr,3), '   Std = ',num2str(Sr,3), '   CoV = ',num2str(CoV,2),  ' Fano = ',num2str(Sr^2/Mr,2) ]);


%%  2 enhancer percent time active model

% The target activation for an individual gene is N/N_tot nuclei
% if driven by a single enhancer this process is at best Poisson, and the
% probability of getting 10% lower or 10% higher than the target is given
% by the incomplete gamma function.

% This premise is wrong. 

gammainc(1,3,'lower'); % chance of getting less than 1 event when on average we get 3 events per window
gammainc(1,2,'lower')^2; 


% When there are two enhancers in the system things become more
% complicated. 
% 

sum(poissrnd(1.1,1,N)<1)/N;


sum(cumsum(-log(rand(1,1000))/7 )<10); % poisson with mean 70

sum(cumsum(-log(rand(1,1000))/4.5 )<10); % poisson with mean 45

clear all; 

N = 100; 
oneEr = .7;
twoEr = oneEr/2;

samples = 1000;
oneE = zeros(1,samples);
twoE = zeros(1,samples);
for n = 1:samples;
    oneE(n) = sum(poissrnd(oneEr,1,N)<1)/N; 
    twoE(n) = sum(poissrnd(twoEr,1,N)+poissrnd(twoEr,1,N)<1)/N;
end
figure(1); clf; subplot(2,1,1); hist(oneE,linspace(0,1,30));
subplot(2,1,2); hist(twoE,linspace(0,1,30)); 

std(oneE)
std(twoE)
