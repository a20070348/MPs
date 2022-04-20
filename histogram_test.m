%data=load('/Users/shuolinxiao/Desktop/vel_scin0160000cz10.txt'); 
%hist(data,1000)
 data2=load('/Users/shuolinxiao/Desktop/SImutaneous_Tracking/vel_scin0180000cz10p2n.dat');
 data3=load('/Users/shuolinxiao/Desktop/SImutaneous_Tracking/vel_scin0180000cz10p3n.dat');
 data4=load('/Users/shuolinxiao/Desktop/SImutaneous_Tracking/vel_scin0180000cz10p4n.dat');
% data5=load('/Users/shuolinxiao/Desktop/SImutaneous_Tracking/vel_scin0160000cz10p5.dat');
 ws4=load('/Users/shuolinxiao/Desktop/SImutaneous_Tracking/vel_scin0180000cz10ws4.dat');
data2p=load('/Users/shuolinxiao/Desktop/SImutaneous_Tracking/vel_scin0180000cz10p2pn.dat');
data4p=load('/Users/shuolinxiao/Desktop/SImutaneous_Tracking/vel_scin0180000cz10p4pn.dat');
h2 = histfit(data2,100,'kernel');
delete(h2(1))
set(h2(2),'color','m')
set(h2(2),'LineStyle','--')
hold on
h3 = histfit(data3,100,'kernel');
delete(h3(1))
set(h3(2),'color','r')
hold on
h4 = histfit(data4,100,'kernel');
set(h4(2),'color','b')
set(h4(2),'LineStyle','-.')
delete(h4(1))
xlim([0 200])
set(gca,'FontSize',24);
title(' Histogram of Fiber Concentration','fontsize',18)
xlabel('$c$','fontsize',36,'Interpreter','latex') 
ylabel('$count$','fontsize',36,'Interpreter','latex')
legend({'horizontally aligned fiber','randomly aligned fiber','fiber settling dependent on local dissipation rate'},'Location','northeast','fontsize',24)
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

