%Sctipt for plotting steady-state values of theta states (theta) while
%changing costs of migration d
%default matching function shape (mtype) is Cobb-Douglas. Can change the
%default by setting optset mtype ='HM'

clc
clear
%optset('match','mtype','HM'); %globally defines matching function type
B=zeros(24,14,'single');
M=zeros(24,1);
i=1;
dm=1;%costs of migration
while dm<2.2
[J,U,W,u,wp,t,p] = country2_ss([0.8,0.8], [1,6.5], [dm,1], [0.2,0.2], [0.3,0.3]);
M(i)=dm;
B(i,:)=[J,U,W,u,wp,t,p];
dm=dm+0.05;
i=i+1;
end
B;
%% 8. plots
theta1=B(:,11);
theta2=B(:,12);
dmig=M;
f=figure(1);
plot(dmig, theta1, 'r');
hold on
plot(dmig, theta2);
hold on
%hx = graph2d.constantline(1.101764195, 'LineStyle','--', 'Color',[.7 .7 .7]);%migration threshhold
%changedependvar(hx,'x');
xlabel('d, costs of migration');
ylabel('theta');
legend('t1','t2')
%legend('unemp1','unemp2','migration threshhold')

hold off