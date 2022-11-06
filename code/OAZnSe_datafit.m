% clear all; 
tau=120e-15;%120fs  脉冲宽度 7-9ns 
ZnSe.wvl=[560,575:25:725,740,750,775,790,800:25:950];%nm 激光能量

L_GaP=0.5e-3;%280um
L_ZnSe=1e-3;%400um

% 
% xlabel('$Z(mm)$','interpreter','latex');
% ylabel('Normalized T');
% % text(40,0.8,'Closed Aperture','FontSize',16);
% xlim([-120 150])



j=0;
for i=ZnSe.wvl
    j=j+1;
    ZnSe.oadata{j}.wavelength=i;
    ZnSe.oadata{j}.data=importdata("../data-202206/ZnSe/"+num2str(i)+"oa.txt").data;
    ZnSe.cadata{j}.wavelength=i;
    ZnSe.cadata{j}.data=importdata("../data-202206/ZnSe/"+num2str(i)+"ca.txt").data;
end
   
%%
for j=1:5
    delta(j)=35*3.375;
    Z{j}=(ZnSe.oadata{j}.data(:,1)-delta(j));%mm
    T{j}=ZnSe.oadata{j}.data(:,4)/mean(ZnSe.oadata{j}.data([1:j*2,(80-j*2):80],4))*0.95;
    T2{j}=ZnSe.cadata{j}.data(:,4)/mean(ZnSe.cadata{j}.data([1:j*2,(80-j*2):80],4))*0.95;
end   
for j=6:7
    delta(j)=35*3.375;
    Z{j}=(ZnSe.oadata{j}.data(:,1)-delta(j));%mm
    T{j}=ZnSe.oadata{j}.data(:,4)/mean(ZnSe.oadata{j}.data([1:j*2,(80-j*2):80],4))*0.9;
    T2{j}=ZnSe.cadata{j}.data(:,4)/mean(ZnSe.cadata{j}.data([1:j*2,(80-j*2):80],4))*0.9;
end   
for j=8:11
    delta(j)=54*2.25;   
    Z{j}=(ZnSe.oadata{j}.data(:,1)-delta(j));%mm
    T{j}=ZnSe.oadata{j}.data(:,4)/mean(ZnSe.oadata{j}.data([1:j*2,(120-j*2):120],4))*0.9;
    T2{j}=ZnSe.cadata{j}.data(:,4)/mean(ZnSe.cadata{j}.data([1:j*2,(120-j*2):120],4))*0.9;
end   
for j=12:16
    delta(j)=22*3.375;
    Z{j}=(ZnSe.oadata{j}.data(:,1)-delta(j));%mm
    T{j}=ZnSe.oadata{j}.data(:,4)/mean(ZnSe.oadata{j}.data([1:j*2,(80-j*2):80],4))*0.8;
    T2{j}=ZnSe.cadata{j}.data(:,4)/mean(ZnSe.oadata{j}.data([1:j*2,(80-j*2):80],4))*0.8;
end   
%%
for j=17:19
    delta(j)=21*3.375;
    Z{j}=(ZnSe.oadata{j}.data(:,1)-delta(j));%mm
    T{j}=ZnSe.oadata{j}.data(:,4)/mean(ZnSe.oadata{j}.data([1:j,(80-j):80],4));
    T2{j}=ZnSe.cadata{j}.data(:,4)/mean(ZnSe.cadata{j}.data([1:j,(80-j):80],4));
%      T3{j}=T2{j}./T{j}
%%
%     subplot(3,1,1)
%     plot(Z{j},T{j},'ok')
%      hold on;
%     %开孔fit 得q0,z0
%     fit=createFit(Z{j},T{j});
%     fitcoefficient=coeffvalues(fit); 
%     q0{j}=2*sqrt(2)*(1-min(T{j}));
%     z0{j}=fitcoefficient(2);%mm
%     r=sqrt(z0{j}*1e-3*(ZnSe.oadata{j}.wavelength)*1e-9/pi); %m   
% %     r=300e-6;%束腰半径27微米 from 张
%     I0=1/(pi*r^2*tau)*1e-6;%焦点处光强J/(m^2s)
% 
%     Toa{j}=1-q0{j}/(2*sqrt(2))*(1./(1+Z{j}.^2/z0{j}^2));%拟合曲线 
%     plot(Z{j},Toa{j},'k')
%     Leff=L_ZnSe;
%     set(gca,'FontSize',16,'FontName','Helvetica','Layer','top');
%     beta(j)=q0{j}/Leff/I0*1e11;%cm/GW
%     hold off;
% 
%     subplot(3,1,2)
%     plot(Z{j},T2{j},'o')
%     set(gca,'FontSize',16,'FontName','Helvetica','Layer','top');
% 
%     subplot(3,1,3)
% %     plot(Z{j},T{j},'o')
%     p{j}=plot(Z{j},T3{j},'ok')
%     hold on;
%     set(gca,'FontSize',16,'FontName','Helvetica','Layer','top');
%     Tpv=max(T3{j})-min(T3{j});
%     fit=CAfit(Z{j}/z0{j},T3{j},Tpv);
%     fitcoa{j}=coeffvalues(fit)   
%     S=fitcoa{j}(1);%q0=beta I0 Leff
%     t=fitcoa{j}(2)
%     n2=(ZnSe.wvl(j)*1e-9/(2*pi))*(Tpv/(0.406*(1-S)^0.27))/(I0*Leff)*1e13%cm2/Gw
%     Tmodel=1+(4*(Tpv/(0.406*(1-S)^0.27)).*(Z{j}/z0{j}-t))./((1+(Z{j}/z0{j}-t).^2).*(9+(Z{j}/z0{j}-t).^2))-0.08
%     plot(Z{j},Tmodel,'k');
%     hold off;
%     saveas(gca,[num2str(ZnSe.oadata{j}.wavelength),'.pdf'])
end  
%% 拟合


for j=1:19

    plot(Z{j},T{j},'ok');
    hold on;       
    fit=createFit(Z{j},T{j});
    fitcoefficient=coeffvalues(fit); 
    q0=2*sqrt(2)*(1-min(T{j}));
%     q0=fitcoefficient(1);%q0=beta I0 Leff
    z0=fitcoefficient(2);%mm
    % q0=2*sqrt(1-min(GaAs.OA{i}));
    r=sqrt(z0*1e-3*(ZnSe.oadata{j}.wavelength)*1e-9/pi); %m   
%     r=300e-6;%束腰半径27微米 from 张
    I0=1/(pi*r^2*tau)*1e-6;%焦点处光强J/(m^2s)  
    T_GaAs=1-q0/(2*sqrt(2))*(1./(1+Z{j}.^2/z0^2));%拟合曲线
    plot(Z{j},T_GaAs);

    Leff=L_ZnSe;
    beta(j)=q0/Leff/I0*1e11;%cm/GW
end
%% 保存数据
en=1240./ZnSe.wvl(1:19)';
head=['energy(eV)','beta(cm/GW)']
data=[en,beta'];
save ZnSe_beta head data
