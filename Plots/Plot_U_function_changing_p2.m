%Sctipt for plotting Values of Unempl state (U) while
%changing productivity in country2 keeping p1 constant
%default matching function shape (mtype) is Cobb-Douglas. Can change the
%default by setting optset mtype ='HM'

clc
clear
optset('match','mtype','CD'); %globally defines matching function type
B=zeros(48,14,'single');
i=1;
prod=1;
while prod<3.4
[J,U,W,u,wp,t,p] = country2_ss_2([0.6,0.6], [1, prod], [1.1,1.1]);
B(i,:)=[J,U,W,u,wp,t,p];
prod=prod+0.05;
i=i+1;
B;
end
%% 8. plots
U_val1=B(:,3);
U_val2=B(:,4);
pr2=B(:,14);
f=figure(1);
plot(pr2, U_val1, 'r');
hold on
plot(pr2, U_val2);
hold on
hx = graph2d.constantline(1.101764195, 'LineStyle','--', 'Color',[.7 .7 .7]);%migration threshhold
changedependvar(hx,'x');
xlabel('p, productivity of country 2');
ylabel('Unemployment state values');
legend('U_value_fun1','U_value_fun_2')
%legend('unemp1','unemp2','migration threshhold')
hold off