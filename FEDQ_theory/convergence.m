%This program is used to calculate the convergence accuracy of the program, that is, to calculate the balance under a series of grid densities
%This program may take a long time
%pre=(0:0.2:1.8)';pre=2.^pre;
pre=(0:0.2:0.4)';pre=2.^pre;

shuju=zeros(length(pre),7);
sss=0;
for k=pre'
    sss=sss+1;
    eq=initial_from_theory(floor(64*k),floor(256*k));
    [eq,error] = FDEQ_theory(eq);
    shuju(sss,:)=error';
    % error(1)=floor(sqrt(eq.geometry.nr*(eq.geometry.t_max-eq.geometry.t_min+1)));%实际网格点密度
    % error(2)=abs(abs_err_psi_Max(end));
    % error(3)=abs(abs_err_psi_mean(end));
    % error(4)=abs(abs_err_psi_var(end));
    % error(5)=err_0(end);
    % error(6)=err_X(end);
    % error(7)=eq.r_range;
end
p2=polyfit(log(shuju(:,1)),log(shuju(:,2)),1);
jingdufit2=exp(polyval(p2,log(shuju(:,1))));
p3=polyfit(log(shuju(:,1)),log(shuju(:,3)),1);
jingdufit3=exp(polyval(p3,log(shuju(:,1))));
p4=polyfit(log(shuju(:,1)),log(shuju(:,4)),1);
jingdufit4=exp(polyval(p4,log(shuju(:,1))));
p7=polyfit(log(shuju(:,1)),log(shuju(:,7)),1);
jingdufit7=exp(polyval(p7,log(shuju(:,1))));


fff=figure;
f2=loglog(shuju(:,1),shuju(:,2),'b*',shuju(:,1),jingdufit2,'b--','LineWidth',2);
hold on
f3=loglog(shuju(:,1),shuju(:,3),'r*',shuju(:,1),jingdufit3,'r--','LineWidth',2);
xlabel('$$\sqrt{n_{r}\times n_{t}}$$','Interpreter','Latex')
ylabel('error of \psi')
legend('Max error',['y \propto x^' '{' num2str(p2(1),2) '}'],'Mean error',['y \propto x^' '{' num2str(p3(1),2) '}'])
fff.Color=[1,1,1];
fff.Units="centimeters";
kuan=fff.Position(3);
gao=fff.Position(4);
kuanfig=6.5;
gaofig=kuanfig/kuan*gao;
set(gcf,'unit','centimeters','Position',[12,2,kuanfig+1.5+0.5,gaofig+1.5+0.5]);  
set(gca,'unit','centimeters','Position',[1.5 1.5 kuanfig gaofig]); 
exportgraphics(gcf,'accuracy.png','Resolution',1000)
