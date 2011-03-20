

%%                      Simplest Enhancer Switching Model

% Alistair Boettiger                                Date Begun: 03/19/11
% Levine Lab            


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
