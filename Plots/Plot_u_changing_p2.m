%Sctipt for plotting steady-state values of unemployment value states (U) while
%changing productivity in country2 keeping p1 constant
%default matching function shape (mtype) is Cobb-Douglas. Can change the
%default by setting optset mtype ='HM'

clc
clear
%optset('match','mtype','HM'); %globally defines matching function type
B=zeros(24,14,'single');
i=1;
prod=1;
while prod<3.4
[J,U,W,u,wp,t,p] = country2_ss([0.8,0.8], [1, prod], [1.1,1.1], [0.2,0.2], [0.2,0.2]);
B(i,:)=[J,U,W,u,wp,t,p];
prod=prod+0.1;
i=i+1;
end
B;
%% 8. plots
unemp1=B(:,7);
unemp2=B(:,8);
pr2=B(:,14);
f=figure(1);
plot(pr2, unemp1, 'r');
hold on
plot(pr2, unemp2);
hold on
%hx = graph2d.constantline(1.101764195, 'LineStyle','--', 'Color',[.7 .7 .7]);%migration threshhold
%changedependvar(hx,'x');
xlabel('p, productivity of country 2');
ylabel('unemployment rates');
legend('unemp1','unemp2')
%legend('unemp1','unemp2','migration threshhold')

hold off


