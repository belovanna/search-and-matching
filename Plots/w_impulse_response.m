%Sctipt for plotting steady-state l.market tightness (theta) while
%changing productivity in country2 keeping p1 constant
%default matching function shape (mtype) is Cobb-Douglas. Can change the
%default by setting optset mtype ='HM'

clc
clear
optset('match','mtype','CD'); %globally defines matching function type
B=zeros(48,14,'single');
i=1;
PR=ones(48,2);
PR(4,2)=3.3;
prod=1;
while i<49
[J,U,W,u,wp,t,p] = country2_ss([0.8,0.8], [1, PR(i,2)], [1.1,1.1],[0.2,0.2],[0.3,0.3]);
B(i,:)=[J,U,W,u,wp,t,PR(i,:)];
i=i+1;
end
%% 8. plots
w1=B(:,9);
w2=B(:,10);
pr2=B(:,14);
f=figure(1);
plot(pr2, w1, 'r');
hold on
plot(pr2, w2);
hold on
%hx = graph2d.constantline(1.139184, 'LineStyle','--', 'Color',[.7 .7 .7]);%migration threshhold
%changedependvar(hx,'x');
xlabel('p, productivity of country 2');
ylabel('wages');
legend('w1','w2')
%legend('unemp1','unemp2','migration threshhold')

hold off