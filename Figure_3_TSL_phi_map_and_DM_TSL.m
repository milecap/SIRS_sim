% Figure_3_TSL_phi_map_and_DM_TSL

clear; clc; close all
fz = 12;
%% --------------------------------------------------------
% Figure 3)a,b,c)
% Contrast map TSL vs Phi for BASL, RESL, and CRESL
% Base file: SLprep_17_MC_part_6.m
% --------------------------------------------------------
% Parameters definition

%Common parameters
gamma = 4258;      %Hz/G
T1 = 1268/1000;    %s
T2 = 100/1000;     %s
T1r = 165/1000;    %s
T2r = 200/1000;    %s
M0 = [0 0 1];      %Initial magnetization

% Tip down / tip up pulse
Tp = 2.5/1000;     %Pulse duration (s)
tds = 10;          %tipdownsamples
dt = Tp/tds;       %time step tip down pulse: defines the time step of the full sequence
time = linspace(0,Tp,tds);
%tdf = zeros(1,tds);     %angular frequency of tip down pulse. (here is static in the simply rotating frame)

alpha = pi/2;
Btd = alpha/(2*pi*gamma*Tp);%amplitude of tip down pulse
tda = 2*pi*gamma*Btd;   %larmor frequency of tip down pulse

e(:,1) = [0,0,1];
e_c = e(:,1);

Vpp = fliplr([146]);    %value in mVpp.
[BNCpp,~] = BNC_calculator(Vpp); %this returns the magnetic field in nT;
aBNC = BNCpp*1E-5;             % Oscillating field amplitude, G
fBNC = 90;
freq = 90;

TSL_v = 70:1:100;
phi_v = -pi:0.1:pi;

for SL_t = 1:3
    fprintf('SL_t = %d \n',SL_t);
    for j = 1:numel(TSL_v)
        TSL = TSL_v(j)/1000;
        sls = ceil(TSL/dt);    % spin lock samples
        
        if SL_t == 1
            sla(1:sls/2) = 2*pi*freq;%Spin lock amplitude in first half
            sla(sls/2+1:sls) = 2*pi*freq;%Spin lock amplitude in second half
        else
            sla(1:sls/2) = 2*pi*freq;%Spin lock amplitude in first half
            sla(sls/2+1:sls) = -2*pi*freq;%Spin lock amplitude in second half
        end
        %slf = zeros(1,sls);     %Spin lock frequency. Static in simply rotating frame
        
        if SL_t == 3
            it = Tp;            %180 Inversion pulse
            ip = 1;             %1=equal to tip down
        else
            it = 0;             %No inversion pulse
            ip = -1;            %1 = oposite to tip down
        end
        
        %B0 inhomogeneity
        dw0 =0;
        if dw0 ~= 0
            sla(1:sls/2) = sqrt(sla(1:sls/2).^2 + dw0^2);
            sla(sls/2+1:sls) = sign(sla(sls/2+1:sls)).*sqrt(sla(sls/2+1:sls).^2 + dw0^2);
        end
        
        for p = 1:length(phi_v)
            phi = phi_v(p);%rand*2*pi; %initial phase of the BNC
            
            %Definition of BNC field
            nca = 2*pi*gamma*aBNC;	%neuronal current field amplitude BNC
            wNC = 2*pi*fBNC;        %frequency of neuronal field BNC
            
            % Call to Bloch simulation
            [~,~,M(1,:)] = Bloch(M0,T1,T2,time,tda,TSL,sla,dw0,it,ip,nca,wNC,phi,T1r,T2r,alpha);
            [~,~,Mi(1,:)] = Bloch(M0,T1,T2,time,tda,TSL,sla,dw0,it,ip,0,0,0,T1r,T2r,alpha);
            
            M = mean(M(:,end),1);
            Mi = mean(Mi(:,end),1);
            
            z(SL_t,j,p) = (M/Mi).*100;
            
            clear M Mi
        end
    end
end

%Final graph
figure(10)
txt = ['BASL ';'RESL ';'CRESL'];
set(gca, 'DefaultFigureRenderer', 'painters');
set(gcf,'color','white');
tlo = tiledlayout(2,3,'TileSpacing','compact');
for SL_t = 1:3
    ht(SL_t) = nexttile; 
    imagesc(TSL_v,phi_v,squeeze(permute(z(SL_t,:,:),[1,3,2])));
    
    grid on;
    [col,num,typ,scheme] = brewermap(20,'GnBu');
    set(ht(SL_t),'Colormap',col,'Clim',[min(z,[],'all') 100])
    
    if SL_t == 1
        yticks([-2,0,2]);
        ylabel('{\phi} (rad)','FontSize',fz)
    else
        set(gca,'yticklabel',[])
    end
    set(gca,'xticklabel',[])

    title(txt(SL_t,:),'FontSize',fz);
end

hcb = colorbar();
colorTitleHandle = get(hcb,'Title');
titleString = '% M/M_0';
set(colorTitleHandle ,'String',titleString);

%% Part 

load('Figure_DM_TSL_prep.mat'); % This file already contains all the necessary files to create Figure 3) b) i,ii,iii)

figure(10)
for i = 1:3
    nexttile;
    grid on; box on; xlim([70 100]);ylim([0.45 1.25]);
    hold on
    plot(TSL,aux(i,:)./100,'o','MarkerSize',4,'MarkerFaceColor',c(1,:),'MarkerEdgeColor',c(1,:));
    plot(TSL,(Y_mean_trig_off(i,:)./Y_mean_off(i,:)),'^','MarkerSize',4,'MarkerFaceColor',c(2,:),'MarkerEdgeColor',c(2,:));
    plot(TSL_v,squeeze(z(i,:,:)),'--','Color','k');
    plot(TSL_v,squeeze(z_p0(i,:,:)),'Color','k');

    xlabel('TSL (ms)','FontSize',fz)

    if i == 1
        ylabel('% M/M_0','FontSize',fz)
        leg = legend({'meas {\phi} = 0', 'meas {\phi} \in [-\pi, \pi]','sim {\phi} = 0','sim {\phi} \in [-\pi, \pi]'},'Location','southwest');
        leg.ItemTokenSize = [17, 16];
        leg.Position;
    else
        set(gca,'yticklabel',[])
    end
end
set(gcf,'color','white')
set(gcf,'Position',[2227         388         734         433])
