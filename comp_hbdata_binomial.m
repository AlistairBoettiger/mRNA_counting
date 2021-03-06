
%% comp_hbdata_binomial

% Alistair Boettiger                    Date Begun: 01/09/11
% Levine Lab                            Last Modified: 01/15/11

%% Updates
% 01/15/11 revised.  can't simply add size of shadow to size of primary for
% effective enhancer size.  The one that measures the gradient has to be
% the one to loop.  
% correlated effects will enter when there is a probability the promoter
% doesn't respond.  


%% Simulate binomial Response across different times. 

figure(4); clf;  

colordef black; set(gcf,'color','k'); 

for cc =11:13; 

if cc == 13;
    time = 18;
    N = 50;
elseif cc == 12;
    time = 5;
    N = 250;
elseif cc == 11;
    time = 5; 
    N =125; 
end

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
 
 col = {'blue','green','red','magenta'}; 
pL1 = [.23,1,.23,1];% probability of S enhancer not forming a loop in shadow alone, primary alone, and Primary + Shadow
pL2 = [1,.07,.07,.07];% " " for primary
a_effs = [.03, .05, .03+.05]/1.9; 
pro_open = 1; 


for e = 1:3% 4; 
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
        
                % Fusion enhancer
        a = a_effs(3); 
          phalf = gammainc(a*cG*D*T,floor(Th),'lower');  % 1 - upper is lower
         [err,thresh] = min( abs(phalf - .5) );
         pint =  1-gammainc(a*C*D*T,floor(thresh),'upper');
        [err,cut] = min( abs(C - low));  % find where lines intersect
        miss_rate3(t) = pint(cut);  
        
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

 
    
        bf(t,:) = binopdf(0:1:N,N,(1-(1-miss_rate1(t))*(1-pL1(e) ))*(1-(1-miss_rate2(t))*(1-pL2(e) ))*pro_open   );
        if e==4
             bf(t,:) = binopdf(0:1:N,N,(1-(1-miss_rate3(t))*(1-pL2(e) ))*pro_open   );   
        end
    end

%     ns = floor(logspace(0,log10(Ts),10));
%     figure(1); clf; plot(bf(ns,:)','r');
%     set(gcf,'color','k');
    
    xbf = linspace(0,1,N+1);
    bfD = sum(bf);
    bfD = bfD/sum(bfD)*N;
    sum(bfD)*1/N
    
    
    figure(4); subplot(3,1,14-cc); hold on; plot(xbf,bfD,'color',col{e},'linewidth',2);
end
legend('no primary','no shadow','control');
title(['cc',num2str(cc),' T = ',num2str(time), 'min  N = ',num2str(N), ' cells']); 
set(gcf,'color','w');
end

%sum(bfD)*1/N


%% Dostatni numbers

V = 100;


%770/(.003*1*4.8)  % time to count all 770 molecules (takes 891 min, error 3%)
% 1/sqrt(X) = .1 -> X = 100 molecules need to count to have 10% error.
% 100/(.003*1*4.8) --> 115 minutes (Thomas numbers)
% 100/(.003 * 4.5 *4.8) --> still 25 minutes, Dostatni numbers
% 100/(.003*10 * 4.5 *4.8) my enhancer size estimate 2.5 minutes


a = .003;
T = 25*60;% 7*60;
D = 4.5;
c = 4.8;
BP =    1/sqrt(D*a*c*T);

    theta = 93;% 25; % 

C = linspace(2, 7, V); 

intp = zeros(1,V);
for k=1:V;
    c= C(k);  % c = 2 
    lambda = a*c*D*T;
    n = 0:100;
    p = lambda.^n./factorial(n)*exp(-lambda);
%      figure(2); clf; 
%      plot(p,'linewidth',4); set(gcf,'color','k'); colordef black;
%      set(gca,'color','k');
%     sum( p(logical(1-isnan(p)))  );


   %  sum(p(1:theta))  % the manual integration method
   %  gammainc(a*c*D*T,floor(theta),'upper') % the cdf method
    intp(k) = 1-gammainc(a*c*D*T,floor(theta),'upper');
end

h1 = figure(3); clf; set(gca,'color','w'); colordef white; 
plot(C,intp); hold on; plot([4.8,4.8],[0,1],'k'); 
plot([4.3,4.3],[0,1],'r'); 
legend('prob detecting > \theta','boundary conc.','10% less than boundary', 'Location','NorthWest' ); 
xlabel('c, molecules / um^3');
ylabel(['probability of seeing >', num2str(theta),' molecules in time T']);
title(['a = ',num2str(a,3), 'um  T = ',num2str(T,3),'s  D = ',num2str(D,2), 'um^2/s', '  1/(DacT)^{1/2} = ' num2str(BP,3)]);
set(gcf,'color','w');



 print('-depsc', '-tiff', [fout,'dostatni_calc_25min.eps']);

%%  multiple binding sites new effective receptor size of primary enhancer

% factor of pi from binding model version of Berg-Purcell

a = .003;
T = 7*60;
c = 4.8;
b = a*24;
m=6;
D = 4.5;

dcM = 1/sqrt(pi*D*c*T).*sqrt(1/(m*a) + 1/(2*b))

%dcD = 1/sqrt(D*c*T*a_eff) = dcM
a_eff = 1/(dcM*sqrt(D*c*T))^2


%%
V = 100;

a = a_eff;
T = 7*60;
D = 4.5;
c = 4.8;

BP =    1/sqrt(D*a*c*T);
C = linspace(2, 7, V); 
intp = zeros(1,V);
for k=1:V;
    c= C(k); 
    lambda = a*c*D*T;
    theta = 430; 
    intp(k) =  1- gammainc(lambda,floor(theta),'upper');
end

h2 = figure(3); clf; plot(C,intp); hold on; plot([4.8,4.8],[0,1],'k'); 
plot([4.3,4.3],[0,1],'r'); 
legend('prob detecting > \theta','boundary conc.','10% less than boundary', 'Location','NorthWest' ); 
xlabel('c, molecules / um^3');
ylabel(['probability of seeing >', num2str(theta),' molecules in time T']);
title(['a = ',num2str(a,3), 'um  T = ',num2str(T,3),'s  D = ',num2str(D,2), 'um^2/s', '  1/(DacT)^{1/2} = ' num2str(BP,3)]);
set(gcf,'color','w');

% saveas(h2,[fout,'aeff_calc.eps'],'eps');

print('-depsc', '-tiff', [fout,'aeff_calc2.eps']);

%% Effective size of shadow

a = .003;
T = 7*60;
c = 4.8;
b = a*20;  % the 3 bcd sites are only 200 bp apart
m=4;
D = 4.5;

dcM = 1/sqrt(pi*D*c*T).*sqrt(1/(m*a) + 1/(2*b))

%dcD = 1/sqrt(D*c*T*a_eff) = dcM
a_eff = 1/(dcM*sqrt(D*c*T))^2


a = a_eff;
T = 7*60;
D = 4.5;
c = 4.8;

BP =    1/sqrt(D*a*c*T);
C = linspace(2, 7, V); 
intp = zeros(1,V);
for k=1:V;
    c= C(k); 
    lambda = a*c*D*T;
    theta = 295; 
    intp(k) =  1- gammainc(lambda,floor(theta),'upper');
end

h2 = figure(3); clf; plot(C,intp); hold on; plot([4.8,4.8],[0,1],'k'); 
plot([4.3,4.3],[0,1],'r'); 
legend('prob detecting > \theta','boundary conc.','10% less than boundary', 'Location','NorthWest' ); 
xlabel('c, molecules / um^3');
ylabel(['probability of seeing >', num2str(theta),' molecules in time T']);
title(['a = ',num2str(a,3), 'um  T = ',num2str(T,3),'s  D = ',num2str(D,2), 'um^2/s', '  1/(DacT)^{1/2} = ' num2str(BP,3)]);
set(gcf,'color','w');

% saveas(h2,[fout,'aeff_calc.eps'],'eps');

print('-depsc', '-tiff', [fout,'shadow_calc.eps']);




a = .003;
T = 7*60;
c = 4.8;
b = a*4;  % the 3 bcd sites are almost back to back
m=3;
D = 4.5;

dcM = 1/sqrt(pi*D*c*T).*sqrt(1/(m*a) + 1/(2*b))

%dcD = 1/sqrt(D*c*T*a_eff) = dcM
a_eff_Driever = 1/(dcM*sqrt(D*c*T))^2


%% Compare to binomoial errors
clear all;

load hbSD_010611

figure(1); clf; 

%     cc11  12  13  14 %
Ncc = [125,250,500,1000]; % expected number of on cells.  
cc = {cc11; cc12; cc13; cc14}; 

for c = 1:4
    mu = zeros(1,G);
    sigma = zeros(1,G);
    bi_sig = zeros(1,G); 
    plot_miss = cell(1,G); 
     for k=1:G;     
         plot_miss{k} = foff{k}(cc{c}{k}) ;%foff{k}(cc14{k});
         mu(k) = mean(plot_miss{k})*Ncc(c);
         sigma(k) = std(plot_miss{k})*Ncc(c);
         bi_sig(k) = mu(k)*( 1- nanmean(plot_miss{k}) );
     end
     rels = logical(1-isnan(sigma));

     corc = corrcoef(bi_sig(rels), sigma(rels)   )  ;
     
     D = nanmedian((sigma - bi_sig)./bi_sig);
     
     disp(sigma);
     disp(bi_sig);
     
     figure(1); subplot(2,2,c); scatter(sigma,bi_sig); 
     lin = [0,bi_sig,1.1*max(bi_sig)];
     hold on; plot(lin,lin,'k');
     xlim([0,1.1*max(sigma)]);
     ylim([0,1.1*max(bi_sig)]);
     title(['hb nucs: ', num2str(Ncc(c)), '  corr= ',num2str(corc(1,2),3), '  dist = ',num2str(D)]); 
 
end

set(gcf,'color','w');


