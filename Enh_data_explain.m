
%%  Figures for hunchback enhancer modeling ms

% 



%%  
N = 200; 
lambda = 20;
f = zeros(1,N+1); y = f; h = f;
for x = 0:N
    f(x+1) = poisspdf(x,lambda);
    y(x+1) = pdf('norm',x,lambda,sqrt(lambda));
    h(x+1) = pdf('bino',x,N,lambda/N);
end

x=0:N; 



    figure(1); clf; plot(log(x), f,'r.');
    hold on; plot(log(x),y,'y.');
   
   figure(2); clf; plot(x,h,'r.'); set(gcf,'color','k');
   set(gca,'fontsize',12);



%% Figure explaining curves from histograms.

% uses data from Plot_DotComps2_hb.m

% data = plot_miss;

x = linspace(0,1,30);  % range and number of bins for histogram
xx = linspace(0,1,pts);
sigma = .1;
G = length(data);
X = cell(1,G); XX = cell(1,G); M = cell(1,G); S = cell(1,G);


for k=1:G
    X{k} = x;
    XX{k} =xx;
    M{k} = method;
    S{k} = sigma;
end



var_hist = cellfun(@nonzeros,data,'UniformOutput',false);
var_hist = cellfun(@hist,var_hist,X,'UniformOutput',false);
N_var = cellfun(@sum,var_hist,'UniformOutput',false);
dist_var = cellfun(@hist2dist,var_hist,X,XX,M,S,'UniformOutput',false);


 set(gcf,'color','k'); colordef black; 

figure(1); clf; hist(data{1},linspace(.0,1,30))
  h = findobj(gca,'Type','patch');
  set(h,'EdgeColor','w','FaceColor','r'); 
hold on; 
y = dist_var{1}/max(dist_var{1})*22;
 plot(xx,y,'color','r','LineWidth',3); hold on;

 
%  figure(2); clf; 

 
 hist(data{3},linspace(.0,1,30))
   h = findobj(gca,'Type','patch');
  alpha(.5);
xlim([-.02,1]);
hold on; 
y = dist_var{3}/max(dist_var{3})*6;
 plot(xx,y,'color','b','LineWidth',3);
 ylim([0,25]);

 
 ylabel('number of embryos');
 xlabel('fraction of missing nuclei'); 
 
 %% Normalized version

figure(1); clf; 
x =linspace(.0,1,30);
d1 = hist(data{1},x);
y = dist_var{1};
dN = max(y); 
d1 = d1/max(d1)*dN; 
bar(x,d1,'EdgeColor','w','FaceColor','r'); 
h = findobj(gca,'Type','patch');
set(h,'EdgeColor','w','FaceColor','r'); 
hold on; 

plot(xx,y,'color','r','LineWidth',3); hold on;

 
 
  

    y = dist_var{3};
    dN = max(y); 
     d2 = hist(data{3},x);
    d2 = d2/max(d2)*dN;
    bar(x,d2);
    h = findobj(gca,'Type','patch');
    alpha(.5);
    xlim([-.02,1]);
  
 plot(xx,y,'color','b','LineWidth',3);
 ylim([0,.5]);

 
 ylabel('probability density');
 xlabel('fraction of missing nuclei'); 