%Sctipt for plotting steady-state values of unemployment value states (u) while
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
while dm<3.4
[J,U,W,u,wp,t,p] = country2_ss([0.8,0.8], [1, 6.3], [dm,dm], [0.2,0.2], [0.3,0.3]);
M(i)=dm;
B(i,:)=[J,U,W,u,wp,t,p];
dm=dm+0.1;
i=i+1;
end
B;
%% 8. plots
unemp1=B(:,7);
unemp2=B(:,8);
dmig=M;
f=figure(1);
plot(dmig, unemp1, 'r');
hold on
plot(dmig, unemp2);
hold on
%hx = graph2d.constantline(1.101764195, 'LineStyle','--', 'Color',[.7 .7 .7]);%migration threshhold
%changedependvar(hx,'x');
xlabel('d, costs of migration');
ylabel('unemployment rates');
legend('unemp1','unemp2')
%legend('unemp1','unemp2','migration threshhold')

hold off