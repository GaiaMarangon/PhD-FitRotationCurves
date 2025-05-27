clear all, close all, clc

%--- PARAMETERS ---------------------------------------------------------

%--- GALAXY
% [nameGal,indexGal] = get_galaxyParams("UGC02953");  % with bulge
[nameGal,indexGal] = get_galaxyParams("NGC5055");   % without bulge
% [nameGal,indexGal] = get_galaxyParams("UGC09037");   % without bulge

%--- MODEL
[rhoNames, factor4pi, nameFactor4pi] = get_modelParams("Exp","4pi");
% [rhoNames, factor4pi, nameFactor4pi] = get_modelParams("TruncatedPlummer","4pi");
% [rhoNames, factor4pi, nameFactor4pi] = get_modelParams("HardBall","4pi");

%--- SPECIFIC COMBO PARAMETERS
QChosen = 1;
nChosen = 4;
xChosen = 37;
mChosen = 0.07;

%--- BEST COMBO PARAMETERS
nrBest = 5;
componentsVisibility = 'on';

%--- LOAD PARAMETERS
load_pathfile = "results/"+nameGal+"/"+rhoNames(1)+"_"+nameFactor4pi+"/";
save_pathfile = "figures/"+nameGal+"_"+rhoNames(1)+"_"+nameFactor4pi+"/";
save_namefile = nameGal+"_"+rhoNames(1)+"_"+nameFactor4pi;

%--- GRAPHICS
colorVec_default = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880];[0.3010 0.7450 0.9330];[0.6350 0.0780 0.1840]];
% figure settings
set(groot,'defaultAxesFontSize',20)     % figures font size
set(groot,'DefaultTextFontSize',20)     % figures font size
set(groot, 'defaultAxesTickLabelInterpreter','latex');  % latex labels
set(groot, 'defaultLegendInterpreter','latex');         % latex labels
set(groot, 'defaultTextInterpreter','latex');           % latex labels

%--- OTHER
percLeft  = 0.2;     %linear fit starts at radius percLeft*rOut
percRight = 0.9;     %linear fit stops at radius percRight*rOut


%--- DATA LOADING -------------------------------------------------------

% load computed combo
load(load_pathfile+"paramTable.mat","paramTable")

%create folders for data saving
createSubfolder(save_pathfile)

% sort combinations by chi2
paramTableSorted = sortrows(paramTable,5);
% adjust number of combinations
nrComb = min(height(paramTable),nrBest);
% use the best combination as a reference
QRef = paramTableSorted(1,1);
nRef = paramTableSorted(1,2);
xRef = paramTableSorted(1,3);
mRef = paramTableSorted(1,4);

% create vectors for parameter variations  
QVec = [];
nVec = [];
xVec = [];
mVec = [];
for i=1:height(paramTableSorted)
    % create vectors of Q
    if (nRef == paramTableSorted(i,2)) && (xRef == paramTableSorted(i,3)) && (mRef == paramTableSorted(i,4))
        QVec = [QVec; paramTableSorted(i,1)]; 
    end
    % create vectors of n
    if (QRef == paramTableSorted(i,1)) && (xRef == paramTableSorted(i,3)) && (mRef == paramTableSorted(i,4))
        nVec = [nVec; paramTableSorted(i,2)]; 
    end
    % create vectors of x
    if (QRef == paramTableSorted(i,1)) && (nRef == paramTableSorted(i,2)) && (mRef == paramTableSorted(i,4))
        xVec = [xVec; paramTableSorted(i,3)]; 
    end
    % create vectors of m 
    if (QRef == paramTableSorted(i,1)) && (nRef == paramTableSorted(i,2)) && (xRef == paramTableSorted(i,3))
        mVec = [mVec; paramTableSorted(i,4)]; 
    end
end
QVec = sort(QVec);
nVec = sort(nVec);
xVec = sort(xVec);
mVec = sort(mVec);



%% 

% %--- SPECIFIC COMBINATION ---------------------------------------------- 
% 
% %load file
% namefile = nameGal+...
%         "_Q"+sprintf('%.0f',QChosen*100)+...
%         "_n"+sprintf('%.0f',nChosen*100)+...
%         "_x"+sprintf('%.0f',xChosen*100)+...
%         "_m"+sprintf('%.0f',mChosen*100)+...
%       "_variables.mat";
% load(load_pathfile+namefile);
% 
% % Plot physical velocities, with all components and (one) best fit
% % NOTE: ups is proportional to mass, so sqrt(ups) is proportional to vel
% figure()
% legend('Location','northeastoutside')
% annotation('textbox',[0.485,0.2, 0.1,0.1],'String',nameGal,'FitBoxToText','on','Interpreter','latex')
% set(gcf,'Position',[100 100 560*1.4 420])
% xlabel("$R$ $(kpc)$")
% ylabel("$V$ $(km/s^2)$")
% hold on
% errorbar(r,vTot,vErr,'Color',[0,0,0,1],'DisplayName',"$V_{tot}^{exp}$")
% plot(r,vTot_pred,'-','Color',colorVec_default(1,:),'LineWidth',1.5,'DisplayName',"$V_{tot}^{pred}$") 
% plot(r,vD,'--','Color',[0,0,0,0.3],'LineWidth',1.5,'DisplayName',"$V_{Disc}$")
% plot(r,vB,'-.','Color',[0,0,0,0.3],'LineWidth',1.5,'DisplayName',"$V_{Bulge}$")
% plot(r,vG,     ':' ,'Color',[0,0,0,0.3],'LineWidth',1.5,'DisplayName',"$V_{Gas}$")
% plot(r,vDM_cmp,'--','Color',colorVec_default(1,:),'LineWidth',1.5,'DisplayName',"$V_{Dark Matter}$")
% % save figure
% saveas(gcf,save_pathfile+"specificCombo_"+save_namefile+".png")
% 
% % Print to screen parameters
% fprintf("Specific fit parameters:\n")
% fprintf("Q = %.2f,\tn = %d,\tx = %.2f,\tm = %.2f\n",QChosen,nChosen,xChosen,mChosen)
% fprintf("chiSquare: %.2f\n",chi2)
% 


%%

%--- BEST COMBINATIONS -------------------------------------------------

% set figure    
figure()
set(gcf,'Position',[100 100 560*1.7 420+40*nrComb])
annotation('textbox',[0.5,0.8, 0.1,0.1],'String',nameGal,'FitBoxToText','on','Interpreter','latex')
xlabel("$R$ $(kpc)$")
ylabel("$V$ $(km/s^2)$")
hold on

% load results for best combos
for i=1:nrComb
    % set combo parameters
    Q = paramTableSorted(i,1);
    n = paramTableSorted(i,2);
    x = paramTableSorted(i,3);
    m = paramTableSorted(i,4);
    %load file
    namefile = nameGal+...
            "_Q"+sprintf('%.0f',Q*100)+...
            "_n"+sprintf('%.0f',n*100)+...
            "_x"+sprintf('%.0f',x*100)+...
            "_m"+sprintf('%.0f',m*100)+...
          "_variables.mat";
    load(load_pathfile+namefile);
    % common settings
    lab1(i) = "$\chi^2$ = "+sprintf('%.2f',chi2)+", Q = "+sprintf('%.2f',Q)+", n = "+sprintf('%.0f',n)+", $\xi$ = "+sprintf('%.0f',x)+", m = "+sprintf('%.2f',m);
    transp = i/nrComb;
    % addo to plot
    p1(i) = plot(r,vTot_pred,'-','Color',[colorVec_default(1,:),transp],'LineWidth',1.5);
    if strcmp(componentsVisibility,"on")
        plot(r,vDM_cmp, '--','Color',[colorVec_default(1,:),transp],'LineWidth',1.5,"HandleVisibility","off")
    end
end
%add experimental data
p(1) = errorbar(r,vTot,vErr,'Color',[0,0,0,1]);     lab(1) = "$V_{Tot}^{Exp}$";
lim = axis;
ylimInf = min(lim(3),0);
%add components
if strcmp(componentsVisibility,"on")
    %add common plots
    p(2) = plot(r,vD,'--','Color',[0,0,0,0.3],'LineWidth',1.5);     lab(2) = "$V_{Disc}$"; j=3;
    if any(vB)
        p(j) = plot(r,vB,'-.','Color',[0,0,0,0.3],'LineWidth',1.5);     lab(j) = "$V_{Bulge}$"; j=j+1;
    end
    p(j) = plot(r,vG,     ':' ,'Color',[0,0,0,0.3],'LineWidth',1.5);    lab(j) = "$V_{Gas}$"; j=j+1;
    p(j) = plot(nan,nan,'-', 'Color','k','LineWidth',1.5);  lab(j) = "$V_{Tot}^{Pred}$"; j=j+1;
    p(j) = plot(nan,nan,'--','Color','k','LineWidth',1.5);  lab(j) = "$V_{DM}$"; 
    lim = axis;
    ylimInf = min(lim(3),0);
    ylim([ylimInf-nrComb*(lim(4)-ylimInf)*0.08,lim(4)*1.05]);
    %add legends
    leg1=legend(p,lab,'Location','northeastoutside');
    ah1=axes('position',get(gca,'position'),'visible','off');
    leg2=legend(ah1,p1,lab1,'Location','south'); 
    set(leg2.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
else
    ylim([ylimInf,lim(4)*1.05]);
    legend([p,p1],[lab,lab1],'Location','south')
end

% save figure
if strcmp(componentsVisibility,"on")
    fileName = save_pathfile+"bestCombo_"+num2str(nrBest)+"_"+save_namefile+".png";
else
    fileName = save_pathfile+"bestCombo_"+num2str(nrBest)+"_off_"+save_namefile+".png";
end
saveas(gcf,fileName)

%%

%--- VARY SINGLE PARAMETERS - VARY Q -------------------------------------

close all

% common setting
textParams = "n = "+string(nRef)+", $\xi$ = "+string(xRef)+", m = "+string(mRef);

% figure setting - eigenfunction
figure(1)
leg=legend('Location','northeast');
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
annotation('textbox',[0.13,0.18, 0.1,0.1],'String',textParams,'FitBoxToText','on','EdgeColor','none','Interpreter','latex')
set(gcf,'Position',[100 100 560*1 420])
xlabel("$r$")
ylabel("$f$")
hold on
% figure setting - potential
figure(2)
leg=legend('Location','southeast');
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
annotation('textbox',[0.17,0.82, 0.1,0.1],'String',textParams,'FitBoxToText','on','EdgeColor','none','Interpreter','latex')
set(gcf,'Position',[100 100 560*1 420])
xlabel("$r$")
ylabel("$\phi$")
hold on
% figure setting - velocity (num)
figure(3)
leg=legend('Location','southeast');
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
annotation('textbox',[0.13,0.82, 0.1,0.1],'String',textParams,'FitBoxToText','on','EdgeColor','none','Interpreter','latex')
set(gcf,'Position',[100 100 560*1 420])
xlabel("$r$")
ylabel("$v$")
hold on
% figure setting - velocity, total (phys)
figure(4)
leg=legend('Location','northeastoutside');
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
annotation('textbox',[0.13,0.82, 0.1,0.1],'String',textParams,'FitBoxToText','on','EdgeColor','none','Interpreter','latex')
set(gcf,'Position',[100 100 560*1.4 420])
xlabel("$R$ $(kpc)$")
ylabel("$V$ $(km/s^2)$")
hold on

% loop through possible Q
for i = 1:length(QVec)
    % fix Q
    Q = QVec(i);
    % load files
    namefile = nameGal+...
        "_Q"   +sprintf('%.0f',Q   *100)+...
        "_n"   +sprintf('%.0f',nRef*100)+...
        "_x"   +sprintf('%.0f',xRef*100)+...
        "_m"   +sprintf('%.0f',mRef*100);
    load(load_pathfile+namefile+"_variables.mat")
    % common setting
    label = "Q = "+string(Q);
    transp = i/(length(QVec));
        % % unit conversion to nonDim: characteristic radius
        % r0D_Num = r0D/r_scale_full;
        % r0B_Num = r0B/r_scale_full;
        % r0G_Num = r0G/r_scale_full;
        % % unit conversion to nonDim: characteristic density
        % aD_Num  = aD*upsD /rho_scale_full;
        % aB_Num  = aB*upsB /rho_scale_full;
        % aG_Num  = aG      /rho_scale_full;
    % plot f_Num(r)
    figure(1)
    plot(rDM_Num,f_Num,'-','Color',[colorVec_default(1,:),transp],'LineWidth',1.5,"DisplayName",label)
        % plot(rDM_Num,aD_Num*exp(-rDM_Num/r0D_Num),'--','Color',[colorVec_default(2,:),transp],'LineWidth',1.5,"DisplayName",label)
        % plot(rDM_Num,aB_Num*exp(-rDM_Num/r0B_Num),'--','Color',[colorVec_default(2,:),transp],'LineWidth',1.5,"DisplayName",label)
        % plot(rDM_Num,aG_Num*exp(-rDM_Num/r0G_Num),'--','Color',[colorVec_default(2,:),transp],'LineWidth',1.5,"DisplayName",label)
    % plot phi_Num(r)
    figure(2)
    plot(rDM_Num,phi_Num,'-','Color',[colorVec_default(1,:),transp],'LineWidth',1.5,"DisplayName",label)
    % plot v_Num(r)
    figure(3)
    plot(rDM_Num,vDM_Num,'-','Color',[colorVec_default(1,:),transp],'LineWidth',1.5,"DisplayName",label)
    % v_Phys(r)
    figure(4)
    plot(r,vTot_pred,'-','Color',[colorVec_default(1,:),transp],'LineWidth',1.5,"DisplayName",label)
end

%add other elemens
figure(4)
errorbar(r,vTot,vErr,'Color',[0,0,0,1],"DisplayName","$V_{Tot}^{Exp}$")
lim = axis;
axis([lim(1:2),-10,lim(4)*1.2])

% save figure - eigenfunction
figure(1)
saveas(gcf,save_pathfile+"varyQ_eigenfunction_"+save_namefile+".png")
% save figure - potential
figure(2)
saveas(gcf,save_pathfile+"varyQ_eigenpotential_"+save_namefile+".png")
% save figure - velocity
figure(3)
saveas(gcf,save_pathfile+"varyQ_eigenvelocity_"+save_namefile+".png")
% save figure - velocity
figure(4)
saveas(gcf,save_pathfile+"varyQ_vTotPhys_"+save_namefile+".png")

%%

%--- VARY SINGLE PARAMETERS - VARY n -------------------------------------

close all

% common setting
textParams = "Q = "+string(QRef)+", $\xi$ = "+string(xRef)+", m = "+string(mRef);

% figure setting - eigenfunction
figure(1)
leg=legend('Location','northeast');
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
annotation('textbox',[0.13,0.18, 0.1,0.1],'String',textParams,'FitBoxToText','on','EdgeColor','none','Interpreter','latex')
set(gcf,'Position',[100 100 560*1 420])
xlabel("$r$")
ylabel("$f$")
hold on
% figure setting - potential
figure(2)
leg=legend('Location','southeast');
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
annotation('textbox',[0.17,0.82, 0.1,0.1],'String',textParams,'FitBoxToText','on','EdgeColor','none','Interpreter','latex')
set(gcf,'Position',[100 100 560*1 420])
xlabel("$r$")
ylabel("$\phi$")
hold on
% figure setting - velocity (num)
figure(3)
leg=legend('Location','southeast');
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
annotation('textbox',[0.13,0.82, 0.1,0.1],'String',textParams,'FitBoxToText','on','EdgeColor','none','Interpreter','latex')
set(gcf,'Position',[100 100 560*1 420])
xlabel("$r$")
ylabel("$v$")
hold on
% figure setting - velocity, total (phys)
figure(4)
leg=legend('Location','northeastoutside');
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
annotation('textbox',[0.13,0.82, 0.1,0.1],'String',textParams,'FitBoxToText','on','EdgeColor','none','Interpreter','latex')
set(gcf,'Position',[100 100 560*1.4 420])
xlabel("$R$ $(kpc)$")
ylabel("$V$ $(km/s^2)$")
hold on

% loop through possible n
for i = 1:length(nVec)
    % fix n
    n = nVec(i);
    % load files
    namefile = nameGal+...
        "_Q"   +sprintf('%.0f',QRef*100)+...
        "_n"   +sprintf('%.0f',n   *100)+...
        "_x"   +sprintf('%.0f',xRef*100)+...
        "_m"   +sprintf('%.0f',mRef*100);
    load(load_pathfile+namefile+"_variables.mat")
    % common setting
    label = "n = "+string(n);
    transp = i/(length(nVec));
    % plot f_Num(r)
    figure(1)
    plot(rDM_Num,f_Num,'-','Color',[colorVec_default(1,:),transp],'LineWidth',1.5,"DisplayName",label)
    % plot phi_Num(r)
    figure(2)
    plot(rDM_Num,phi_Num,'-','Color',[colorVec_default(1,:),transp],'LineWidth',1.5,"DisplayName",label)
    % plot v_Num(r)
    figure(3)
    plot(rDM_Num,vDM_Num,'-','Color',[colorVec_default(1,:),transp],'LineWidth',1.5,"DisplayName",label)
    % v_Phys(r)
    figure(4)
    plot(r,vTot_pred,'-','Color',[colorVec_default(1,:),transp],'LineWidth',1.5,"DisplayName",label)
end

%add other elemens
figure(1)
lim = axis;
axis([lim(1:2),lim(3)*1.8,lim(4)])
figure(2)
lim = axis;
axis([lim(1:2),lim(3),-lim(3)*0.1])
figure(4)
errorbar(r,vTot,vErr,'Color',[0,0,0,1],"DisplayName","$V_{Tot}^{Exp}$")
lim = axis;
axis([lim(1:2),-10,lim(4)*1.2])

% save figure - eigenfunction
figure(1)
saveas(gcf,save_pathfile+"varyn_eigenfunction_"+save_namefile+".png")
% save figure - potential
figure(2)
saveas(gcf,save_pathfile+"varyn_eigenpotential_"+save_namefile+".png")
% save figure - velocity
figure(3)
saveas(gcf,save_pathfile+"varyn_eigenvelocity_"+save_namefile+".png")
% save figure - velocity
figure(4)
saveas(gcf,save_pathfile+"varyn_vTotPhys_"+save_namefile+".png")

%%

%--- VARY SINGLE PARAMETERS - VARY x -------------------------------------

close all

%-- VARY x --------------------------------------------------------------

% common setting
textParams = "Q = "+string(QRef)+", n = "+string(nRef)+", m = "+string(mRef);

% figure setting - eigenfunction
figure(1)
leg=legend('Location','northeast');
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
annotation('textbox',[0.13,0.18, 0.1,0.1],'String',textParams,'FitBoxToText','on','EdgeColor','none','Interpreter','latex')
set(gcf,'Position',[100 100 560*1 420])
xlabel("$r$")
ylabel("$f$")
hold on
% figure setting - potential
figure(2)
leg=legend('Location','southeast');
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
annotation('textbox',[0.17,0.82, 0.1,0.1],'String',textParams,'FitBoxToText','on','EdgeColor','none','Interpreter','latex')
set(gcf,'Position',[100 100 560*1 420])
xlabel("$r$")
ylabel("$\phi$")
hold on
% figure setting - velocity (num)
figure(3)
leg=legend('Location','southeast');
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
annotation('textbox',[0.13,0.82, 0.1,0.1],'String',textParams,'FitBoxToText','on','EdgeColor','none','Interpreter','latex')
set(gcf,'Position',[100 100 560*1 420])
xlabel("$r$")
ylabel("$v$")
hold on
% figure setting - velocity, total (phys)
figure(4)
leg=legend('Location','northeastoutside');
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
annotation('textbox',[0.13,0.82, 0.1,0.1],'String',textParams,'FitBoxToText','on','EdgeColor','none','Interpreter','latex')
set(gcf,'Position',[100 100 560*1.4 420])
xlabel("$R$ $(kpc)$")
ylabel("$V$ $(km/s^2)$")
hold on

% loop through possible x
for i = 1:length(xVec)
    % fix x
    x = xVec(i);
    % load files
    namefile = nameGal+...
        "_Q"   +sprintf('%.0f',QRef*100)+...
        "_n"   +sprintf('%.0f',nRef*100)+...
        "_x"   +sprintf('%.0f',x   *100)+...
        "_m"   +sprintf('%.0f',mRef*100);
    load(load_pathfile+namefile+"_variables.mat")
    % common setting
    label = "$\xi$ = "+string(x);
    transp = i/(length(xVec));
    % plot f_Num(r)
    figure(1)
    plot(rDM_Num,f_Num,'-','Color',[colorVec_default(1,:),transp],'LineWidth',1.5,"DisplayName",label)
    % plot phi_Num(r)
    figure(2)
    plot(rDM_Num,phi_Num,'-','Color',[colorVec_default(1,:),transp],'LineWidth',1.5,"DisplayName",label)
    % plot v_Num(r)
    figure(3)
    plot(rDM_Num,vDM_Num,'-','Color',[colorVec_default(1,:),transp],'LineWidth',1.5,"DisplayName",label)
    % v_Phys(r)
    figure(4)
    plot(r,vTot_pred,'-','Color',[colorVec_default(1,:),transp],'LineWidth',1.5,"DisplayName",label)
end

%add other elemens
figure(4)
errorbar(r,vTot,vErr,'Color',[0,0,0,1],"DisplayName","$V_{Tot}^{Exp}$")
lim = axis;
axis([lim(1:2),-10,lim(4)*1.2])

% save figure - eigenfunction
figure(1)
saveas(gcf,save_pathfile+"varyx_eigenfunction_"+save_namefile+".png")
% save figure - potential
figure(2)
saveas(gcf,save_pathfile+"varyx_eigenpotential_"+save_namefile+".png")
% save figure - velocity
figure(3)
saveas(gcf,save_pathfile+"varyx_eigenvelocity_"+save_namefile+".png")
% save figure - velocity
figure(4)
saveas(gcf,save_pathfile+"varyx_vTotPhys_"+save_namefile+".png")



%%

%--- VARY SINGLE PARAMETERS - VARY m -------------------------------------

close all

% common setting
textParams = "Q = "+string(QRef)+", n = "+string(nRef)+", $\xi$ = "+string(xRef);

% figure setting - eigenfunction
figure(1)
leg=legend('Location','northeast');
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
annotation('textbox',[0.13,0.18, 0.1,0.1],'String',textParams,'FitBoxToText','on','EdgeColor','none','Interpreter','latex')
set(gcf,'Position',[100 100 560*1 420])
xlabel("$r$")
ylabel("$f$")
hold on
% figure setting - potential
figure(2)
leg=legend('Location','southeast');
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
annotation('textbox',[0.17,0.82, 0.1,0.1],'String',textParams,'FitBoxToText','on','EdgeColor','none','Interpreter','latex')
set(gcf,'Position',[100 100 560*1 420])
xlabel("$r$")
ylabel("$\phi$")
hold on
% figure setting - velocity (num)
figure(3)
leg=legend('Location','southeast');
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
annotation('textbox',[0.13,0.82, 0.1,0.1],'String',textParams,'FitBoxToText','on','EdgeColor','none','Interpreter','latex')
set(gcf,'Position',[100 100 560*1 420])
xlabel("$r$")
ylabel("$v$")
hold on
% figure setting - velocity, total (phys)
figure(4)
leg=legend('Location','northeastoutside');
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
annotation('textbox',[0.13,0.82, 0.1,0.1],'String',textParams,'FitBoxToText','on','EdgeColor','none','Interpreter','latex')
set(gcf,'Position',[100 100 560*1.4 420])
xlabel("$R$ $(kpc)$")
ylabel("$V$ $(km/s^2)$")
hold on

% loop through possible m
for i = 1:length(mVec)
    % fix m
    m = mVec(i);
    % load files
    namefile = nameGal+...
        "_Q"   +sprintf('%.0f',QRef*100)+...
        "_n"   +sprintf('%.0f',nRef*100)+...
        "_x"   +sprintf('%.0f',xRef*100)+...
        "_m"   +sprintf('%.0f',   m*100);
    load(load_pathfile+namefile+"_variables.mat")
    % common setting
    label = "m = "+string(m);
    transp = i/(length(mVec));
    % plot f_Num(r)
    figure(1)
    plot(rDM_Num,f_Num,'-','Color',[colorVec_default(1,:),transp],'LineWidth',1.5,"DisplayName",label)
    % plot phi_Num(r)
    figure(2)
    plot(rDM_Num,phi_Num,'-','Color',[colorVec_default(1,:),transp],'LineWidth',1.5,"DisplayName",label)
    % plot v_Num(r)
    figure(3)
    plot(rDM_Num,vDM_Num,'-','Color',[colorVec_default(1,:),transp],'LineWidth',1.5,"DisplayName",label)
    % v_Phys(r)
    figure(4)
    plot(r,vTot_pred,'-','Color',[colorVec_default(1,:),transp],'LineWidth',1.5,"DisplayName",label)
end
%add other elemens
figure(4)
errorbar(r,vTot,vErr,'Color',[0,0,0,1],"DisplayName","$V_{Tot}^{Exp}$")
lim = axis;
axis([lim(1:2),-10,lim(4)*1.2])

% save figure - eigenfunction
figure(1)
saveas(gcf,save_pathfile+"varym_eigenfunction_"+save_namefile+".png")
% save figure - potential
figure(2)
saveas(gcf,save_pathfile+"varym_eigenpotential_"+save_namefile+".png")
% save figure - velocity
figure(3)
saveas(gcf,save_pathfile+"varym_eigenvelocity_"+save_namefile+".png")
% save figure - velocity
figure(4)
saveas(gcf,save_pathfile+"varym_vTotPhys_"+save_namefile+".png")

%%

%-- PLOT CHI2 SURFACE PLOTS --------------------------------------------
figure()
set(gcf,'Position',[100 100 560*2 420*2])
% Q,n
subplot(3,3,1)
[X,Y] = meshgrid(QVec,nVec);
Z = zeros(size(X));
for i=1:size(X,1)
    for j=1:size(X,2)
        idxRow = find(ismembertol(paramTable(:,1:4), [round(X(i,j),2),round(Y(i,j)),xRef,mRef],0.001,'ByRows',true));
        if size(idxRow,1)==0
            Z(i,j) = NaN;
        else
            Z(i,j) = paramTable(idxRow,5);
        end
    end
end
%remove NaNs
X(:,sum(isnan(Z),1)>1) = [];
Y(:,sum(isnan(Z),1)>1) = [];
Z(:,sum(isnan(Z),1)>1) = [];
X(sum(isnan(Z),2)>1,:) = [];
Y(sum(isnan(Z),2)>1,:) = [];
Z(sum(isnan(Z),2)>1,:) = [];
%plot
surf(X,Y,Z)
xlabel('$Q$')
ylabel('$n$')
zlabel('$\chi^2$')

% Q,x
subplot(3,3,4)
[X,Y] = meshgrid(QVec,xVec);
Z = zeros(size(X));
for i=1:size(X,1)
    for j=1:size(X,2)
        idxRow = find(ismembertol(round(paramTable(:,1:4),3), [round(X(i,j),2),nRef,round(Y(i,j)),mRef],0.001,'ByRows',true));
        if size(idxRow,1)==0
            Z(i,j) = NaN;
        else
            Z(i,j) = paramTable(idxRow,5);
        end
    end
end
%remove NaNs
X(:,sum(isnan(Z),1)>1) = [];
Y(:,sum(isnan(Z),1)>1) = [];
Z(:,sum(isnan(Z),1)>1) = [];
X(sum(isnan(Z),2)>1,:) = [];
Y(sum(isnan(Z),2)>1,:) = [];
Z(sum(isnan(Z),2)>1,:) = [];
%plot
surf(X,Y,Z)
xlabel('$Q$')
ylabel('$\xi$')
zlabel('$\chi^2$')

% Q,m
subplot(3,3,7)
[X,Y] = meshgrid(QVec,mVec);
Z = zeros(size(X));
for i=1:size(X,1)
    for j=1:size(X,2)
        idxRow = find(ismembertol(round(paramTable(:,1:4),3), [round(X(i,j),2),nRef,xRef,round(Y(i,j),2)],0.001,'ByRows',true));
        if size(idxRow,1)==0
            Z(i,j) = NaN;
        else
            Z(i,j) = paramTable(idxRow,5);
        end
    end
end
%remove NaNs
X(:,sum(isnan(Z),1)>1) = [];
Y(:,sum(isnan(Z),1)>1) = [];
Z(:,sum(isnan(Z),1)>1) = [];
X(sum(isnan(Z),2)>1,:) = [];
Y(sum(isnan(Z),2)>1,:) = [];
Z(sum(isnan(Z),2)>1,:) = [];
%plot
surf(X,Y,Z)
xlabel('$Q$')
ylabel('$m$')
zlabel('$\chi^2$')

% n,x
subplot(3,3,5)
[X,Y] = meshgrid(nVec,xVec);
Z = zeros(size(X));
for i=1:size(X,1)
    for j=1:size(X,2)
        idxRow = find(ismembertol(round(paramTable(:,1:4),3), [QRef,round(X(i,j)),round(Y(i,j)),mRef],0.001,'ByRows',true));
        if size(idxRow,1)==0
            Z(i,j) = NaN;
        else
            Z(i,j) = paramTable(idxRow,5);
        end
    end
end
%remove NaNs
X(:,sum(isnan(Z),1)>1) = [];
Y(:,sum(isnan(Z),1)>1) = [];
Z(:,sum(isnan(Z),1)>1) = [];
X(sum(isnan(Z),2)>1,:) = [];
Y(sum(isnan(Z),2)>1,:) = [];
Z(sum(isnan(Z),2)>1,:) = [];
%plot
surf(X,Y,Z)
xlabel('$n$')
ylabel('$\xi$')
zlabel('$\chi^2$')

% n,m
subplot(3,3,8)
[X,Y] = meshgrid(nVec,mVec);
Z = zeros(size(X));
for i=1:size(X,1)
    for j=1:size(X,2)
        idxRow = find(ismembertol(round(paramTable(:,1:4),3), [QRef,round(X(i,j)),xRef,round(Y(i,j),2)],0.001,'ByRows',true));
        if size(idxRow,1)==0
            Z(i,j) = NaN;
        else
            Z(i,j) = paramTable(idxRow,5);
        end
    end
end
%remove NaNs
X(:,sum(isnan(Z),1)>1) = [];
Y(:,sum(isnan(Z),1)>1) = [];
Z(:,sum(isnan(Z),1)>1) = [];
X(sum(isnan(Z),2)>1,:) = [];
Y(sum(isnan(Z),2)>1,:) = [];
Z(sum(isnan(Z),2)>1,:) = [];
%plot
surf(X,Y,Z)
xlabel('$n$')
ylabel('$m$')
zlabel('$\chi^2$')

% x,m
subplot(3,3,9)
[X,Y] = meshgrid(xVec,mVec);
Z = zeros(size(X));
for i=1:size(X,1)
    for j=1:size(X,2)
        idxRow = find(ismembertol(round(paramTable(:,1:4),3), [QRef,nRef,round(X(i,j)),round(Y(i,j),2)],0.001,'ByRows',true));
        if size(idxRow,1)==0
            Z(i,j) = NaN;
        else
            Z(i,j) = paramTable(idxRow,5);
        end
    end
end
%remove NaNs
X(:,sum(isnan(Z),1)>1) = [];
Y(:,sum(isnan(Z),1)>1) = [];
Z(:,sum(isnan(Z),1)>1) = [];
X(sum(isnan(Z),2)>1,:) = [];
Y(sum(isnan(Z),2)>1,:) = [];
Z(sum(isnan(Z),2)>1,:) = [];
%plot
surf(X,Y,Z)
xlabel('$\xi$')
ylabel('$m$')
zlabel('$\chi^2$')

saveas(gcf,save_pathfile+"chi2_"+save_namefile+".png")

%%

%--- PLOT REFERENCE EXPERIMENTAL FEATURES -------------------------------

%load file (best combo)
namefile = nameGal+...
        "_Q"+sprintf('%.0f',QRef*100)+...
        "_n"+sprintf('%.0f',nRef*100)+...
        "_x"+sprintf('%.0f',xRef*100)+...
        "_m"+sprintf('%.0f',mRef*100)+...
      "_variables.mat";
load(load_pathfile+namefile);

%compute DM experim (ks/s)
vDMsqr_exp = vTot.*abs(vTot) - (vG.*abs(vG) + vB.*abs(vB) + vD.*abs(vD) );
vDM_exp = sqrt(abs(vDMsqr_exp)) .* sign(vDMsqr_exp);
%retrieve reference data: 
[~,idxStart] = min(abs( r-percLeft *r(end) ));
[~,idxStop]  = min(abs( r-percRight*r(end) ));
P = polyfit(r(idxStart:idxStop),vDM_exp(idxStart:idxStop),1);
slope = P(1);
interc = P(2);
for i=3:length(r) %not from beginnning to avoid spurious effects
    if vDM_exp(i) > slope*r(i)+interc
        Rin = r(i);
        Vin = Rin*slope+interc;
        break
    end
end
Rout = r(end);
Vout = vDM_exp(end);
%plot of reference points of vDM_exp
figure()
xlabel('$R$ $(kpc)$')
ylabel('$V$ $(km/s)$')
leg=legend('Location','southeast');
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.8]));
textParams = nameGal + ", Q = "+string(QRef);
annotation('textbox',[0.17,0.82, 0.1,0.1],'String',textParams,'FitBoxToText','on','EdgeColor','none','Interpreter','latex')
set(gcf,'Position',[100 100 560 420])
hold on
plot(r,vDM_exp,'--+k','DisplayName',"$V_{DM}^{Exp}$")
plot(r(idxStart:idxStop),slope*r(idxStart:idxStop)+interc,"-",'Color',colorVec_default(1,:),'LineWidth',1.5,"DisplayName","V = "+string(round(slope,2))+"R+"+ string(round(interc,2)));
plot(Rin,Vin,'x','Color',colorVec_default(2,:),'LineWidth',1.5,'MarkerSize',10,'DisplayName',"$(\tilde{R}_{0},\tilde{V}_{0})$")
plot(Rout,Vout,'x','Color',colorVec_default(3,:),'LineWidth',1.5,'MarkerSize',10,'DisplayName',"$(R_{last},V_{last})$")
lim = axis;
axis([lim(1:2),lim(3),lim(4)*1.1])

saveas(gcf,save_pathfile+"refElemsVDMExp_"+save_namefile+".png")




