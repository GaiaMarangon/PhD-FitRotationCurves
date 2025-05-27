clear all, close all, clc

%--- DESCRIPTION --------------------------------------------------------
% This code tests a uniform grid of parameters, without estimates.
% Given reference parameters (QRef,nRef,xRef,mRef), it varies parameters one
% at a time as specified by (QVec,nVec,xVec,mVec) and computes the
% corresponding combination. 

%--- PARAMETERS ---------------------------------------------------------

%--- GALAXY
% [nameGal,indexGal] = get_galaxyParams("UGC02953");  % with bulge
[nameGal,indexGal] = get_galaxyParams("NGC5055");   % without bulge
% [nameGal,indexGal] = get_galaxyParams("UGC09037");   % without bulge    

%--- MODEL
[rhoNames, factor4pi, nameFactor4pi] = get_modelParams("Exp","4pi");
% [rhoNames, factor4pi, nameFactor4pi] = get_modelParams("TruncatedPlummer","4pi");
% [rhoNames, factor4pi, nameFactor4pi] = get_modelParams("HardBall","4pi");

%--- PARAMETERS TO BE TESTED
% set depending on previous results
% Default reference values:
QRefDef = 0.8;
nRefDef = 7;
xRefDef = 35;
mRefDef = 0.12;

%--- CONSTANTS ----------------------------------------------------------

%--- FIDUCIAL VALUES
upsDFid = 0.5;     % (Ms/Ls) 
upsBFid = 0.7;     % (Ms/Ls)

%--- SOLVER PARAMETERS (TOLERANCES)
epsPhi = 1e-015;        %for regularizing 1/r at r=0
tolInt = 1e-08;         %for finding const eigval, in internal iter
tolExt = 1e-10;         %for finding const eigval, in external iter
maxiterInt = 100;       %maximum number of internal iterations
maxiterExt = 100;       %maximum number of external iterations

%--- PHYSICAL CONSTANTS 
C_lp = 1.62;
C_mp = 2.18;
C_pc = 3.09;
C_ms = 1.99;
C_g  = 6.67;
C_c  = 3.00;
% units conversion for physical constants
lp_kpc = C_lp / C_pc;
mp_9ms = C_mp*10/C_ms;
mp_56g = C_mp;
G_astr = C_g*C_ms/C_pc*1000;

%--- PATH FOR DATA SAVING
pathfile = "results/"+nameGal+"/"+rhoNames(1)+"_"+nameFactor4pi+"/";
paramTableFile = "paramTable.mat";

%--- PATH FOR DATA LOADING
inputFileRV  = "input/Rotmod_LTG/"+nameGal+"_rotmod.dat";
inputFileRho = "input/densityParams.dat";




%--- DATA LOADING -------------------------------------------------------
%load experimental velocity curve (v in [km/s], r in [kpc])
data = importdata(inputFileRV).data;
r = data(:,1);
vTot = data(:,2);
vErr = data(:,3);
vD_lum = data(:,5);
vB_lum = data(:,6);
vG_lum = data(:,4);

%load density parameters (a in [10^9 Ls/kpc^3], r0 in [kpc])
densParams = importdata(inputFileRho).data;
r0D = densParams(indexGal,1);
aD  = densParams(indexGal,2);
r0B = densParams(indexGal,3);
aB  = densParams(indexGal,4);
r0G = densParams(indexGal,5);
aG  = densParams(indexGal,6);

%load available parameter table (if any)
if isfile(pathfile+paramTableFile) 
    load(pathfile+paramTableFile)
    if isempty(paramTable)
        fprintf("Warning: paramTable empty, using default values\n\n")
        warningParam = true;
    else
        warningParam = false;
    end
else 
    fprintf("Warning: paramTable not found, using default values\n\n")
    warningParam = True;
    paramTable = [];
end
%set reference values
if warningParam
    QRef = QRefDef;
    nRef = nRefDef;
    xRef = xRefDef;
    mRef = mRefDef;
else
    % sort combinations by chi2
    paramTableSorted = sortrows(paramTable,5);
    % set reference parameters
    QRef = paramTableSorted(1,1);
    nRef = paramTableSorted(1,2);
    xRef = paramTableSorted(1,3);
    mRef = paramTableSorted(1,4);
end

%set parameter variations
QVec = [QRef-0.1, QRef-0.05, QRef, QRef+0.05, QRef+0.1]; 
nVec = [nRef-2, nRef-1, nRef, nRef+1, nRef+2];   
xVec = [xRef-2, xRef-1, xRef, xRef+1, xRef+2];
mVec = [mRef-0.02, mRef-0.01, mRef, mRef+0.01, mRef+0.02];

% %--- SET GRID, OPTION ONE ----------------------------------------------
% % parameters are varied one at a time
% paramGrid = zeros(length(QVec)+length(nVec)+length(xVec)+length(mVec), 4);
% for i=1:length(QVec)
%     paramGrid(i,:) = [QVec(i),nRef,xRef,mRef];
% end
% for i=1:length(nVec)
%     paramGrid(length(QVec)+i,:) = [QRef,nVec(i),xRef,mRef];
% end
% for i=1:length(xVec)
%     paramGrid(length(QVec)+length(nVec)+i,:) = [QRef,nRef,xVec(i),mRef];
% end
% for i=1:length(mVec)
%     paramGrid(length(QVec)+length(nVec)+length(xVec)+i,:) = [QRef,nRef,xRef,mVec(i)];
% end

%--- SET GRID, OPTION TWO ----------------------------------------------
% all combinations of parameters are computed
paramGrid = zeros(length(QVec)*length(nVec)*length(xVec)*length(mVec), 4);
for iQ=1:length(QVec)
    for in=1:length(nVec)
        for ix=1:length(xVec)
            for im=1:length(mVec)
                paramGrid(length(nVec)*length(xVec)*length(mVec)*(iQ-1)+length(xVec)*length(mVec)*(in-1)+length(mVec)*(ix-1)+im,:) = [QVec(iQ),nVec(in),xVec(ix),mVec(im)];
            end
        end
    end
end

%create folders for data saving
createSubfolder(pathfile)


%--- FITTING PROCEDURE --------------------------------------------------

for i=1:length(paramGrid)

    % set parameters
    Q = paramGrid(i,1);
    n = paramGrid(i,2);
    x = paramGrid(i,3);
    m = paramGrid(i,4);

    % check if combination is already computed
    isComputed = ismembertol(paramTable(:,1:4),paramGrid(i,:),0.001,'ByRows',true);
    if any(isComputed) 
        % combination is alreay computed
        fprintf("Computed combo: Q = "+string(Q) +", n = "+string(n) +", x = "+string(x) +", m = "+string(m)+"\n")
        fprintf("Computed %d-th combo out of %d\n\n",i,height(paramGrid));
    else
        % compute results for the combination

        % namefile for saving
        namefile = nameGal+...
            "_Q" +sprintf('%.0f',Q*100)+...
            "_n" +sprintf('%.0f',n*100)+...
            "_x" +sprintf('%.0f',x*100)+...
            "_m" +sprintf('%.0f',m*100);
    
        %retrieve scales
        upsD = Q*upsDFid;
        upsB = Q*upsBFid;
        %compute velocity components (km/s)
        vD = sqrt(upsD)*vD_lum;
        vB = sqrt(upsB)*vB_lum;
        vG =            vG_lum;
        %compute total Mass (10^9 Ms)
        ML = 8*pi*(upsD*aD*r0D^3 + upsB*aB*r0B^3 + aG*r0G^3);
        %compute scales
        r_scale   = lp_kpc*mp_9ms*(mp_56g)^2/ML /factor4pi;
        rho_scale = ML^4 /mp_9ms^3 /mp_56g^6 /lp_kpc^3 *factor4pi^3;
        v_scale   = sqrt(4*pi)*C_c*100*ML/mp_56g/mp_9ms *sqrt(factor4pi);
        %scale factor, numerical to physical units
        r_scale_full   = r_scale /(m^2*x);
        rho_scale_full = rho_scale *m^6*x^4;
        v_scale_full   =v_scale *m*x;
       
        %fix discretization step
        expDom95 = 10*n^2+17*n+8; 
        if n==0 
            h = expDom95*0.0012;
        elseif n==1
            h = expDom95*0.0008;
        elseif n==2
            h = expDom95*0.0007;
        elseif n==3
            h = expDom95*0.0006;
        elseif n==4
            h = expDom95*0.0006;
        else   %(better not beyond n=10)
            h = expDom95*0.00055;
        end            
                            
        % unit conversion to nonDim: characteristic radius
        r0D_Num = r0D/r_scale_full;
        r0B_Num = r0B/r_scale_full;
        r0G_Num = r0G/r_scale_full;
        % unit conversion to nonDim: characteristic density
        aD_Num  = aD*upsD /rho_scale_full;
        aB_Num  = aB*upsB /rho_scale_full;
        aG_Num  = aG      /rho_scale_full;
        % convert from normOne to 4piNorm notation for solver
        r0D_4pi = r0D_Num / (4*pi);
        r0B_4pi = r0B_Num / (4*pi);
        r0G_4pi = r0G_Num / (4*pi);
        aD_4pi  = aD_Num  * (4*pi)^4;
        aB_4pi  = aB_Num  * (4*pi)^4;
        aG_4pi  = aG_Num  * (4*pi)^4;
        %opening file for solver log 
        fileID = fopen(pathfile+namefile+"_nthSolver.txt",'w');
        % set parameters for phi(r) external contribution
        phiParams_4pi = [r0D_4pi,aD_4pi; r0B_4pi,aB_4pi; r0G_4pi,aG_4pi];
        %solve problem
        tic
        [eigval4pi,r4pi,phi4pi,f4pi] = nthSolver(n,h,phiParams_4pi,rhoNames,epsPhi,tolInt,tolExt,maxiterInt,maxiterExt,fileID);
        elapsedTime = toc;
        %print to file computation time
        fprintf('\nElapsed time: %f\n',elapsedTime);
        fprintf(fileID,'\nElapsed time: %f\n',elapsedTime);
        fclose(fileID);
        %convert 4piNorm notation to NormOne notation
        [eigval_Num, f_Num, phi_Num, rDM_Num, ~]=Norm4piToNormOne_extSource(eigval4pi, f4pi, phi4pi, r4pi, []);
        
        % compute DM velocity (numerical units)
        [~,vDM_Num] = massVelCurv(rDM_Num,f_Num,epsPhi,n,'n','n');
        % convert to physical units
        rDM = rDM_Num*r_scale_full; 
        vDM = vDM_Num*v_scale_full;
        % compute Dark Matter velocity at tabulated radial positions, via interpolation 
        vDM_cmp = interp1(rDM,vDM,r,'linear','extrap');
        % TEMPORARY: plot vDM predicted, in physical units
        figure(1)
        hold on
        plot(r,vDM_cmp,'DisplayName',"Q = "+string(Q) +", n = "+string(n) +", x = "+string(x) +", m = "+string(m))
        legend('Location','northeastoutside')
        set(gcf,'Position',[100 100 560*1.7 420])
        % combine with experimental velocity components to get the total predicted velocity
        vTotSqr_pred = vG.*abs(vG) + vB.*abs(vB) + vD.*abs(vD) + vDM_cmp.*abs(vDM_cmp);
        vTot_pred = sqrt(abs(vTotSqr_pred)) .* sign(vTotSqr_pred);
        
        % compute chi2 metric
        chi2 = sum( (vTot-vTot_pred).^2./(vErr.^2) ) / (length(vTot)-4);
        % add accepted combination with chi2 to table
        paramTable = [paramTable; [paramGrid(i,:), chi2]];
        
        % save results to .mat file
        save(pathfile+namefile+"_variables.mat",...
            "eigval_Num","rDM_Num","f_Num","phi_Num","vDM_Num", ...
            "Q","upsD","upsB","n","x","m","ML",...
            "r0D","r0B","r0G","aD","aB","aG",...
            "r0D_Num","r0B_Num","r0G_Num","aD_Num","aB_Num","aG_Num",...
            "r_scale_full","v_scale_full","rho_scale_full",...
            "rDM","vDM", ...
            "r","vErr","vD","vB","vG","vDM_cmp","vTot","vTot_pred","chi2");
        % print progress to command window
        fprintf("Computed combo: Q = "+string(Q) +", n = "+string(n) +", x = "+string(x) +", m = "+string(m)+"\n")
        fprintf("Computed %d-th combo out of %d\n\n",i,height(paramGrid));
        % update list of computed combinations (create a different one if paramTable was initially empty)
        if warningParam
            save(pathfile+paramTableFile(1:end-4)+"_uniform.mat",'paramTable')
        else
            save(pathfile+paramTableFile,'paramTable')
        end
    end
end
