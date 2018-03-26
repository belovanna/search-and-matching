%Sctipt for plotting steady-state values of unemployment value states (u) while
%changing productivity in country2 keeping p1 constant
%default matching function shape (mtype) is Cobb-Douglas. Can change the
%default by setting optset mtype ='HM'

clc
clear
optset('match','mtype','CD'); %globally defines matching function type
B=zeros(24,7,'single');
i=1;
prod=1;
while prod<3.4
[J,U,W,u,wp,t,p] = country_ss(0.8, prod, 0.2, 0.2);
B(i,:)=[J,U,W,u,wp,t,p];
prod=prod+0.1;
i=i+1;
end
B;
%% 8. plots
unemp=B(:,4);
pr2=B(:,7);
f=figure(1);
plot(pr2, unemp, 'r');

%hx = graph2d.constantline(1.101764195, 'LineStyle','--', 'Color',[.7 .7 .7]);%migration threshhold
%changedependvar(hx,'x');
xlabel('p, productivity of country 2');
ylabel('unemployment rates');
legend('unemp changing p')
%legend('unemp1','unemp2','migration threshhold')

hold off