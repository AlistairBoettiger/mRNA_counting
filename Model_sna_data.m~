

%%                          Model_Kr_kni_data.m

% Alistair Boettiger                    Date Begun: 01/17/11
% Levine Lab                            Last Modified: 01/19/11
%                                           

%% Simulate Kr and kni data
% Adapted from simulation of hb data written in comp_hbdata_binomial.m
% 
%% Updates
%
%% sna data


% inputs

 names = {'no primary 22','no shadow 22','control 22','no primary 30','no shadow 30','control 30'};
 col = {'blue','green','red','cyan','yellow','magenta'}; 

 Ls = .03; % probability of shadow fail looping
 Lp = .01; % prob primary
 pL1 = [Ls,1,Ls,4*Ls,1,4*Ls];% probability of S enhancer not forming a loop in shadow alone, primary alone, and Primary + Shadow
pL2 = [1,Lp,Lp,1,4*Lp,4*Lp];% " " for primary
a_effs = [.02, .02]; 
pro_open = .9; 



figure(4); clf;  

colordef black; set(gcf,'color','k'); 

cc =14; 
time = 30; % time in minutes
N = 20; 
 
V = 200; % number of concentration points to check
stp = 20; 
Ts = time*60/stp; % number of time points to check


D = 4.5;% diffusion rate of bcd (according to Dostatni)
cG = 4.8; % c Gregor

C = linspace(2,7,V);  % range of concentrations to explore; 
Th = 1:1000; % range of thresholds number of molecules to explore
      
low = cG-cG*.1; % 10% less than threshold bcd concentration
high = cG + cG*.1; % 10% more than threshold bcd concentration

miss_rate1 = zeros(1,Ts); % Primary, diffusion limited miss rate
miss_rate2 = zeros(1,Ts);  % shadow, diffusion limited miss rate
miss_rate3 = zeros(1,Ts);  % shadow, diffusion limited miss rate
 



for e = 1:6 
    bf = zeros(Ts,N+1);
    for t = 1:Ts
        T = t*stp;% T = 1*60;
        
        % Primary enhancer
        a = a_effs(1); 
        phalf = gammainc(a*cG*D*T,floor(Th),'lower');  % 1 - upper is lower
         [err,thresh] = min( abs(phalf - .5) );
         pint =  1-gammainc(a*C*D*T,floor(thresh),'upper');
        [err,cut] = min( abs(C - low));  % find where lines intersect
        miss_rate1(t) = pint(cut);       
        
        % Shadow enhancer
        a = a_effs(2); 
          phalf = gammainc(a*cG*D*T,floor(Th),'lower');  % 1 - upper is lower
         [err,thresh] = min( abs(phalf - .5) );
         pint =  1-gammainc(a*C*D*T,floor(thresh),'upper');
        [err,cut] = min( abs(C - low));  % find where lines intersect
        miss_rate2(t) = pint(cut);    
        
        
%       % KEEP THIS.  USEFUL CODE!        
%         figure(1); clf; plot(C,pint,'linewidth',2); hold on; 
%         plot([low,low],[0,1],'c');
%         hold on; plot([high,high],[0,1],'r');
%         figure(1); hold on; plot([C(cut)], [pint(cut)],'r.')
%         
%         vert = [C(1:cut)',pint(1:cut)'];  vert = [vert; C(cut),0;0,0];
%         fac = [1:cut+2];
%         
%         patch('Faces',fac,'Vertices',vert,'Facecolor',[0,.5,1]);
%         alpha(.7);
        
 %   end

%     figure(2); clf; plot(miss_rate, 'k.'); 
%     bvar = N*miss_rate.*(1-miss_rate);
%     figure(2); clf; plot(bvar); 
%     figure(3); clf; hist(miss_rate,0:.01:1); xlim([0,1]);

    
         bf(t,:) = binopdf(0:1:N,N,1-(1- (1- (1-miss_rate1(t))*(1-pL1(e)) )*(1- (1-miss_rate2(t))*(1-pL2(e)) )     )*pro_open  );
    end

%     ns = floor(logspace(0,log10(Ts),10));
%     figure(1); clf; plot(bf(ns,:)','r');
%     set(gcf,'color','k');
    
    xbf = linspace(0,1,N+1);
    bfD = sum(bf);
    bfD = bfD/sum(bfD)*N;
    mean(bfD)
    
    
    figure(4); subplot(2,1,ceil(e/3.1)); hold on; plot(xbf,bfD,'color',col{e},'linewidth',2);
    legend(names);
end

title(['cc',num2str(cc),' T = ',num2str(time), 'min  N = ',num2str(N), ' cells']); 



%sum(bfD)*1/N