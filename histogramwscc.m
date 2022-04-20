clear;
clc;

ws4scc=load('/Users/shuolinxiao/Desktop/SImutaneous_Tracking/fiber_with_small_size/vel_scin0180000cz10ws4scc.dat');
% ws4n=load('/Users/shuolinxiao/Desktop/SImutaneous_Tracking/fiber_with_small_size/vel_scin0180000cz10ws4n.dat');
h2 = histfit(ws4scc,100,'kernel');
delete(h2(1))
set(h2(2),'color','m')
set(h2(2),'LineStyle','--')
hold on
% h3 = histfit(ws4n,100,'kernel');
% delete(h3(1))
% set(h3(2),'color','b')
% set(h3(2),'LineStyle','-')
%xlim([0 0.07])
set(gca,'FontSize',24);
title(' Histogram of Fiber Settling Velocity','fontsize',18)
xlabel('$w_{s}(m/s)$','fontsize',36,'Interpreter','latex') 
ylabel('$count$','fontsize',36,'Interpreter','latex')
legend({'fiber settling dependent on local dissipation rate in CBL'},'Location','northeast','fontsize',24)
% h5 = hist(data5,100);
% hist(ws4,100)
% h6 = hist(data2p,100);
% hold on
% h7 = hist(data4p,100);
%  plot(h6)
%  hold on 
%  plot(h7)
% plot(h2)
% hold on 
%  plot(h3)
% hold on 
%  plot(h4,'black')
 %hold on 
% plot(h5)

