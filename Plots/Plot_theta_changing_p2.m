%Sctipt for plotting steady-state l.market tightness (theta) while
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
[J,U,W,u,wp,t,p] = country2_ss([0.8,0.8], [1, prod], [1.1,1.1],[0.2,0.2],[0.3,0.3]);
B(i,:)=[J,U,W,u,wp,t,p];
prod=prod+0.05;
i=i+1;
end
%% 8. plots
theta1=B(:,11);
theta2=B(:,12);
pr2=B(:,14);
f=figure(1);
plot(pr2, theta1, 'r');
hold on
plot(pr2, theta2);
hold on
xlabel('p, productivity of country 2');
ylabel('tightness');
legend('theta1','theta2')
%legend('unemp1','unemp2','migration threshhold')

hold off