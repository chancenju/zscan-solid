clear all; 

tau=120e-15;%8ns  脉冲宽度 7-9ns 
GaP.wvl=[575:25:725,740,750,775,790,800:25:950];%nm 激光能量
L_GaP=0.5e-3;%
L_ZnSe=1e-3;%400um
% 引入GaP Z-scan 数据
j=0;
for i=GaP.wvl
    j=j+1;
    GaP.oadata{j}.wavelength=i;
    GaP.oadata{j}.data=importdata("..\data-202206\GaP\"+num2str(i)+"oa.txt").data;
    GaP.cadata{j}.wavelength=i;
    GaP.cadata{j}.data=importdata("..\data-202206\GaP\"+num2str(i)+"ca.txt").data;
end

%%
delta(1:4)=58*2.25;
delta(5:7)=59*2.25;
delta([8:10])=56*2.25;
delta([11:15])=52*2.25;
delta([16:18])=52*2.25;
for j=1:18
    Z{j}=GaP.oadata{j}.data(:,1)-delta(j);%mm
    T{j}=GaP.oadata{j}.data(:,4)/mean(GaP.oadata{j}.data([1:j*2,(120-j*2):120],4))*0.99;
    T2{j}=GaP.cadata{j}.data(:,4)/mean(GaP.cadata{j}.data([1:j*2,(120-j*2):120],4))*0.88;
    normalT{j}=T{j}/mean(T{j}(1:10));
    plot(Z{j},T{j},'ok');
    hold on;       
    fit=createFit(Z{j},T{j});
    fitcoefficient=coeffvalues(fit); 
    q0=2*sqrt(2)*(1-min(T{j}));
%     q0=fitcoefficient(1);%q0=beta I0 Leff
    z0=fitcoefficient(2);%mm
    % q0=2*sqrt(1-min(GaAs.OA{i}));
    r=sqrt(z0*1e-3*(GaP.oadata{j}.wavelength)*1e-9/pi); %m   
%     r=300e-6;%束腰半径27微米 from 张
    I0=1/(pi*r^2*tau)*1e-6;%焦点处光强J/(m^2s)  
    T_GaAs=1-q0/(2*sqrt(2))*(1./(1+Z{j}.^2/z0^2));%拟合曲线
    plot(Z{j},T_GaAs);

    Leff=L_GaP;
    beta(j)=q0/Leff/I0*1e11;%cm/GW
%%

end
%%
en=1240./GaP.wvl(1:18)';
head=['energy(eV)','beta(cm/GW)']
data=[en,beta'];
save GaP_beta head data
%%


semilogy(en,beta','ok','MarkerFaceColor','k','linewidth',2)
hold on;
x=en/2.8;
js=30*(2*x'-1).^(3/2)./(2*x').^5;
semilogy(en,js,'--k','LineWidth',2)
set(gcf,'Position',[400,100  ,1000,600]);
set(gca,'FontSize',16,'FontName','Helvetica','Layer','top');
xlabel('photon energy(eV)');
ylabel('\beta_{TPA}(cm/GW)');
ylim([0.1 25])
xlim([1,2])
xticks([1,1.25,1.5,1.75,2])
set(gcf,'Units','inches');
pos=get(gcf,'Position');
set(gcf,"PaperPositionMode","auto","PaperUnits","inches","PaperSize",[pos(3),pos(4)])
print(gcf,'OAGaP.pdf','-dpdf','-r0')
close(gcf)