% clear all
tau=120e-15;%8ns  脉冲宽度 7-9ns 
GaP.wvl=[575:25:725,740,750,775,790,800:25:950];%nm 激光能量
% GaAs.wvl={'1200','1250','1300','1350','1400','1450','1500','1550','1600','1600'}
% GaAs.alpha={10^0.772,10^0.827,10^0.793,10^0.811,10^0.832};%cm^-1 
% GaAs.alpha=10^0.830707;%cm^-1 
L_GaP=0.5e-3;%280um
L_ZnSe=1e-3;%400um
% for i=GaP.wvl
%     1
% end
% CASi{i}={['Closed Aperture',num2str(i),'.txt']};
% OASi{i}={['Open Aperture',num2str(i),'.txt']};
% Si.CASi{i}=importdata(string(CASi{i})).data;
% Si.OASi{i}=importdata(string(OASi{i})).data;
% z=Si.CASi{1}(:,1)-52;
% set(gcf,'position',[0 40 600 800])
% 
% Si.OA{i}=Si.OASi{i}(:,4)/Si.OASi{i}(1,4);
% Si.fit{i}=createFit(z,Si.OA{i})
% fitcoefficient{i}=coeffvalues(Si.fit{i});
% q0=fitcoefficient{i}(1);%q0=beta I0 Leff 美女
% z0=fitcoefficient{i}(2);  
% % q0=2*sqrt(1-min(Si.OA{i}));  
% r=sqrt(z0*1e-3*wvl/pi);    
% % r=27e-6;%束腰半径27微米 from 张
% I0=str2num(Si.E0{i})/(pi*r^2*tau)*1e-6;%焦点处光强J/(m^2s)
% T_Si=1-q0/(2*sqrt(2))*(1./(1+z.^2/z0^2));%拟合曲线
% Leff=(1-exp(-Si.alpha*L_Si))/Si.alpha;%低吸收近似
% beta{i}=q0/Leff/I0;%m/W
% subplot(3,1,1)
% plot(z,Si.OA{i},z,T_Si,'linewidth',2);
% % beta=(1-0.17)*2*sqrt(2)/((40e-6/(3.14*27e-6*27e-6*7e-9))*280e-6)*1e11;
% 
% 
% 
% % plot(z,Si.OA{i},'linewidth',2)
% legend('Experiment','Theory','FontSize',10);
% xlabel('$Z(mm)$','interpreter','latex');
% ylabel('Normalized T');
% text(20,0.7,'Open Aperture','FontSize',16);
% betatext=['beta=',num2str(beta{i}),'m/W']
% % text(-60,0.9,betatext,'FontSize',16);
% set(gca,'FontSize',16,'FontName','Helvetica','Layer','top');
% hold on;
% subplot(3,1,2)
% Si.CA{i}=Si.CASi{i}(:,4)/(Si.CASi{i}(1,4));
% plot(z,Si.CA{i},'linewidth',2)
% xlabel('$Z(mm)$','interpreter','latex');
% ylabel('Normalized T');
% text(20,0.85,'Closed Aperture','FontSize',16);
% set(gca,'FontSize',16,'FontName','Helvetica','Layer','top');
% hold on;
% subplot(3,1,3)
% Si.COA{i}=Si.CASi{i}(:,5)/(Si.CASi{i}(1,5));
% DeltaTpv=max(Si.COA{i})-min(Si.COA{i});
% S=0.4;
% DeltaPhi=DeltaTpv/(0.406*(1-S)^0.27);
% TCOA=1-4*(z/z0)*(DeltaPhi-0.1)./((1+z.^2/z0^2).*(9+z.^2/z0^2));%拟合曲线
% n2{i}=-wvl*DeltaPhi/(2*pi*I0*Leff);%m^2/W
% n2text=['n2=',num2str(n2{i}),'m^2/W']
% 
% 
% 
% 
% plot(z+5,Si.COA{i},'linewidth',2)
% xlabel('$Z(mm)$','interpreter','latex')
% ylabel('Normalized T')
% text(30,1.3,'CA/OA','FontSize',16)
% % text(-60,0.2,n2text,'FontSize',16);
% set(gca,'FontSize',16,'FontName','Helvetica','Layer','top');
% hold on;
% title=['Z-scan of Silicon@',Si.E0{i},'uJ']
% suptitle(title)
% end
j=0;
for i=GaP.wvl
    j=j+1;
    GaP.oadata{j}.wavelength=i;
    GaP.oadata{j}.data=importdata("../data-202206/GaP/"+num2str(i)+"oa.txt").data;
    GaP.cadata{j}.wavelength=i;
    GaP.cadata{j}.data=importdata("../data-202206/GaP/"+num2str(i)+"ca.txt").data;
end

%%
for j=1:4
    delta(j)=58*2.25;
end   
for j=5:7
    delta(j)=59*2.25;
end   
for j=8:10
    delta{j}=56*2.25;
end   
for j=11:15
    delta{j}=52*2.25;
end   
for j=16:18
    delta{j}=52*2.25;
end  
for j=1:18
    Z=(GaP.oadata{j}.data(:,1)-delta{j});%mm
    T=GaP.oadata{j}.data(:,4)/mean(GaP.oadata{j}.data([1:j*2,(120-j*2):120],4));
    plot(Z,T,'o')
    hold on;
    fit=createFit(Z,T);
    fitcoefficient=coeffvalues(fit); 
    q0=2*sqrt(2)*(1-min(T));
%     q0=fitcoefficient(1);%q0=beta I0 Leff
    z0=fitcoefficient(2);%mm
    % q0=2*sqrt(1-min(GaAs.OA{i}));
    r=sqrt(z0*1e-3*(GaP.oadata{j}.wavelength)*1e-9/pi); %m   
%     r=300e-6;%束腰半径27微米 from 张
    I0=1/(pi*r^2*tau)*1e-6;%焦点处光强J/(m^2s)
    T_GaAs=1-q0/(2*sqrt(2))*(1./(1+Z.^2/z0^2));%拟合曲线
    plot(Z,T_GaAs);

    Leff=L_GaP;
    beta(j)=q0/Leff/I0*1e11;%cm/GW
        % hold on;
    % z=GaAs.CAGaAs{1}(:,1)-125;
    % set(gcf,'position',[0 40 600 800])
    % 
    % GaAs.OA{i}=GaAs.OAGaAs{i}(:,4)/mean(GaAs.OAGaAs{i}(1:20,4));
    % GaAs.fit{i}=createFit(z,GaAs.OA{i})
    % fitcoefficient{i}=coeffvalues(GaAs.fit{i});
    % q0=fitcoefficient{i}(1);%q0=beta I0 Leff
    % z0=fitcoefficient{i}(2);
    % % q0=2*sqrt(1-min(GaAs.OA{i}));
    % r=sqrt(z0*1e-3*str2num(GaAs.wvl{i})*1e-9/pi);    
    % % r=27e-6;%束腰半径27微米 from 张
    % I0=300/(pi*r^2*tau)*1e-6;%焦点处光强J/(m^2s)
    % T_GaAs=1-q0/(2*sqrt(2))*(1./(1+z.^2/z0^2));%拟合曲线
    % % Leff=(1-exp(-GaAs.alpha{i}*L_GaAs))/GaAs.alpha{i};%低吸收近似
end
% %%
% for j=1:18
%     Z=(GaP.oadata{j}.data(:,1)-delta{j});%mm
%     T=GaP.oadata{j}.data(:,4)/mean(GaP.oadata{j}.data([1:j*2,(120-j*2):120],4));
%     plot(Z,T,'o')
%     hold on;
%     q0=2*sqrt(2)*(1-min(T));
%     fit=createFit(Z,T);
%     fitcoefficient=coeffvalues(fit); 
% %     q0=fitcoefficient(1);%q0=beta I0 Leff
%     z0=fitcoefficient(2);%mm
%     % q0=2*sqrt(1-min(GaAs.OA{i}));
%     r=sqrt(z0*1e-3*(GaP.oadata{j}.wavelength)*1e-9/pi); %m   
% %     r=300e-6;%束腰半径27微米 from 张
%     I0=1/(pi*r^2*tau)*1e-6;%焦点处光强J/(m^2s)
%     T_GaAs=1-q0/(2*sqrt(2))*(1./(1+Z.^2/z0^2));%拟合曲线
%     plot(Z,T_GaAs);
% 
%     Leff=L_GaP;
%     beta(j)=q0/Leff/I0*1e11;%cm/GW
% end
% plot(GaP.wvl(1:18),beta)
% b(i)=beta{i}
% subplot(3,1,1)
% plot(z,GaAs.OA{i},'b.','linewidth',2);
% 
% % plot(z,GaAs.OA{i},'b.',z,T_GaAs,'linewidth',2);
% % beta=(1-0.17)*2*sqrt(2)/((40e-6/(3.14*27e-6*27e-6*7e-9))*280e-6)*1e11;
% 
% 
% 
% % plot(z,Si.OA{i},'linewidth',2)
% legend('Experiment','Theory','FontSize',10);
% xlabel('$Z(mm)$','interpreter','latex');
% ylabel('Normalized T');
% % text(40,0.9,'Open Aperture','FontSize',16);
% betatext=['beta=',num2str(beta{i}),'m/W']
% xlim([-125 125])
% % text(-120,0.9,betatext,'FontSize',16);
% set(gca,'FontSize',16,'FontName','Helvetica','Layer','top');
% hold on;
% 
% 
% subplot(3,1,2)
% GaAs.CA{i}=GaAs.CAGaAs{i}(:,4)/mean(GaAs.CAGaAs{i}(1:20,4));
% plot(z,GaAs.CA{i},'b.','linewidth',2)
% xlabel('$Z(mm)$','interpreter','latex');
% ylabel('Normalized T');
% % text(40,0.8,'Closed Aperture','FontSize',16);
% xlim([-125 125])
% set(gca,'FontSize',16,'FontName','Helvetica','Layer','top');
% hold on;
% 
% 
% subplot(3,1,3)
% GaAs.COA{i}=GaAs.CAGaAs{i}(:,5)/mean(GaAs.CAGaAs{i}(1:20,5));
% GaAs.COAfit{i}=CAfit(z/z0,GaAs.COA{i});
% fitcoa{i}=coeffvalues(GaAs.COAfit{i})
% DeltaPhi3=fitcoa{i}(1);%q0=beta I0 Leff
% % z0=fitcoefficient{i}(2)
% DeltaPhi5=fitcoa{i}(2);
% % DeltaPhi5=0.5;
% DeltaTpv=-max(GaAs.COA{i})+min(GaAs.COA{i});
% S=0.25;
% DeltaPhi2=DeltaTpv/(0.406*(1-S)^0.27);
% TCOA=1+4*(z/z0)*(DeltaPhi3)./((1+z.^2/z0^2).*(9+z.^2/z0^2))+4*2*(z/z0)*DeltaPhi5./((1+z.^2/z0^2).^2.*(25+z.^2/z0^2));%拟合曲线
% TCOA2=1+4*(z/z0)*(DeltaPhi2)./((1+z.^2/z0^2).*(9+z.^2/z0^2));
% n2{i}=wvl*DeltaPhi3/(2*pi*I0*Leff);%m^2/W
% n(i)=n2{i};
% n2p{i}=wvl*DeltaPhi2/(2*pi*I0*Leff);%m^2/W
% np(i)=n2p{i};
% n2text=['n2_{3+5}=',num2str(n2{i}),'m^2/W'];
% n2textp=['n2_3=',num2str(n2p{i}),'m^2/W'];
% 
% plot(z,GaAs.COA{i},'b.','linewidth',2)
% % plot(z,TCOA,'k',z,TCOA2,z,GaAs.COA{i},'b.','linewidth',2)
% legend('3rd','3rd+5th')
% xlabel('$Z(mm)$','interpreter','latex');
% ylabel('Normalized T');
% % text(60,0.9,'CA/OA','FontSize',16);
% % text(-120,0.8,n2text,'FontSize',14);
% % text(-120,0.9,n2textp,'FontSize',14);
% xlim([-125 125])
% ylim([0.7 1.4])
% set(gca,'FontSize',16,'FontName','Helvetica','Layer','top');
% hold on;
% title=['Z-scan of Si@',GaAs.wvl{i},'nm'];
% suptitle(title)
% end
% i=3
% subplot(3,1,1)
% GaAs.OA{i}=GaAs.OAGaAs{i}(:,4)/max(GaAs.OAGaAs{i}(:,4));
% plot(z,GaAs.OA{i})
% subplot(3,1,2)
% GaAs.CA{i}=GaAs.CAGaAs{i}(:,4)/(GaAs.CAGaAs{i}(1,4));
% plot(z,GaAs.CA{i})
% subplot(3,1,3)
% GaAs.COA{i}=GaAs.CAGaAs{i}(:,5)/(GaAs.CAGaAs{i}(1,5));
% plot(z,GaAs.COA{i})
% zscan.i.GaAs=data.GaAs(:,4)/data.GaAs(1,4);%Normalized transmission
% z=zscan.z;%距离
% i=zscan.i.Si;%透射率

% L=L_Si;
% fit=createFit(z,i)
% %beta I
% fitcoefficient=coeffvalues(fit);
% q0=fitcoefficient(1);%q0=beta I0 Leff
% z0=fitcoefficient(2);  
% % r=sqrt(z0*1e-3*wvl/pi);
% r=27e-6;%束腰半径27微米 from 张
% I0=E0/(pi*r^2*tau);%焦点处光强
% T_GaAs=1-q0/(2*sqrt(2))*(1./(1+z.^2/z0^2));%拟合曲线
% Leff=(1-exp(-alpha*L))/alpha;;%低吸收近似
% beta=q0/Leff/I0*1e11%m/W推到cm/GW
% plot(z,i,'o',z,T_GaAs,'k','Linewidth',2);
% % beta=(1-0.17)*2*sqrt(2)/((40e-6/(3.14*27e-6*27e-6*7e-9))*280e-6)*1e11;
% beta=['$\beta$=',num2str(beta),'cm/GW']; 
% xlabel('$Z(mm)$','interpreter','latex')
% ylabel('Normalized T')
% text(6,0.3,beta,'interpreter','latex','FontSize',16)
% title('TPA of Silicon')
% set(gca,'FontSize',16,'FontName','Helvetica','Layer','top');
% beta={'3.73','2.86','2.44','1.47','0','0'};
% E0=[4.9,10.2,24.6,50]
% n2=[1.79,1.22,1.15,0.5]
% beta2=[3.73,2.86,2.44,1.47]
%%
% for i=t
%     plot(i,length(GaP.oadata{1,i}.data(:,1)),'o')
% hold on;
% end
% plot(wavelenth,np,'k-o','MarkerSize',4,'LineWidth',3,'MarkerFaceColor','k')
% % legend('3rd','3rd+5th')
% xlabel('Wavelength (nm)');
% text(40,5,'tau=8ns','FontSize',16)
% ylabel('n_2(m^2/W)');
% % ylim([0 6])
% % ylabel('\beta(m/W)');
% text(-60,0.9,betatext,'FontSize',16);
% set(gca,'FontSize',16,'FontName','Helvetica','Layer','top');
% title('n_2 of InP')
