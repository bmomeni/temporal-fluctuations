clear

infile = 'pHCheckC_fpH0.2_EnsembleRobustCoexistenceCS_RndExpSpt0LVsa_BC_Ngen100_fp20_Ne10000_rndseed7239';
load(strcat(infile,'.mat'))

crec = 0; % counter for recording different cases
for ne = 1:Ne
    if (DCS(ne)<0.1)&&(Ncxst(ne)>1)
        tic
%         disp(ne)
        crec = crec+1;
        N = NS(ne);
        Nc = Ncxst(ne);
        indx = 1:N;
        SSindx = indx(CompS(1:N,ne)>1e-6);
        pHr(crec) = pH(ne);
        
        cip = ciSS(SSindx,SSindx,ne);
        cntSpecies(crec,1:Ns0) = hist(Cxst(1:Ncxst(ne),ne),1:Ns0);
        cntIntType(crec,1) = sum(sum(cip<0))-Nc; % number of negative interaction coefficients
        cntIntType(crec,2) = sum(sum(cip>0)); % number of positive interaction coefficients
        CompDev(crec) = max(DCd(ne,:));
        CompDevAll(crec,1:Nd) = DCd(ne,:);
        pHdCase(crec) = max(pHdd(ne,:));
        pHdCaseAll(crec,1:Nd) = pHdd(ne,:);
    end
end

FracPosInt = cntIntType(:,2)./sum(cntIntType,2);
figure
plot(pHdCase(FracPosInt>0),CompDev(FracPosInt>0),'.')
xlabel('pH fluctuation amplitude')
ylabel('Composition deviation')

selCoop = (FracPosInt == 0.5);
selComp = (FracPosInt == 0);
DevComp_Comp = CompDev(selComp);
DevComp_Coop = CompDev(selCoop);
DevCompColl = [DevComp_Comp DevComp_Coop];
GDevComp = [ones(1,length(DevComp_Comp)) 2*ones(1,length(DevComp_Coop))];
figure
boxplot(DevCompColl,GDevComp,'width',0.2)
ylabel('Composition deviation')
p = ranksum(DevComp_Comp,DevComp_Coop);
text(1.5,0.8,num2str(p))
text(1.1,0.12,num2str(mean(DevComp_Comp)))
text(1.7,0.06,num2str(mean(DevComp_Coop)))

NpHdbin = 10;
pHdbin = linspace(min(pHdCase),max(pHdCase),NpHdbin+1);
pHdCaseAlllin = reshape(pHdCaseAll,1,crec*Nd);
CompDevAlllin = reshape(CompDevAll,1,crec*Nd);
Gbin = zeros(1,crec*Nd);
for cpb = 1:NpHdbin
    selbin = (pHdCaseAlllin>pHdbin(cpb))&(pHdCaseAlllin<=pHdbin(cpb+1));
    DevOverpHm(cpb) = mean(CompDevAlllin(selbin));
    DevOverpHsd(cpb) = std(CompDevAlllin(selbin));
    Gbin(selbin) = cpb;
end
figure
errorbar(pHdbin(1:NpHdbin)+(pHdbin(2)-pHdbin(1))/2,DevOverpHm,DevOverpHsd)
xlabel('pH fluctuation amplitude')
ylabel('Composition deviation')
legend('Robust communities','Sensitive communities')

figure
boxplot(CompDevAlllin,Gbin,'positions',[0 pHdbin(1:NpHdbin)+(pHdbin(2)-pHdbin(1))/2-0.05],'outliersize',0.1,'width',0.02)
xlabel('pH fluctuation amplitude')
ylabel('Composition deviation')

for cb = 1:crec
    RNDindx = ceil(Ns0*rand(1,Nsp));
    cir = ci(RNDindx,RNDindx);
    cntIntType0(cb,1) = sum(sum(cir<0))-Nsp; % number of negative interaction coefficients
    cntIntType0(cb,2) = sum(sum(cir>0)); % number of positive interaction coefficients
end
FracPosInt0 = cntIntType0(:,2)./sum(cntIntType0,2);
figure
plot(FracPosInt0,FracPosInt,'.')
hold on
plot([0 0.5],[0 0.5],'r')
xlabel('Initial fraction of facilitation')
ylabel('Final fraction of facilitation')
xlim([0 0.5])
ylim([0 0.5])

figure
bar(sum(cntSpecies))
xlabel('Species')
ylabel('Persistance Frequency')

NpHbin = 10;
pHbin = linspace(min(pHr),max(pHr),NpHbin+1);
for cpb = 1:NpHbin
    selbin = (pHr>pHbin(cpb))&(pHr<=pHbin(cpb+1));
    SpeciesFreq(:,cpb) = sum(cntSpecies(selbin,:));
end
figure
semilogy(pHbin(1:NpHbin)+(pHbin(2)-pHbin(1))/2,1/crec*SpeciesFreq')
xlabel('pH')
ylabel('Frequency of persistance')
legend('1821','1828','1850','1867','1989','1839')

%% Sensitive
infile = 'pHCheck_fpH0.2_SensitiveEnsembleRobustCoexistenceCS_RndExpSpt0LVsa_BC_Ngen100_fp20_Ne10000_rndseed7239';
load(strcat(infile,'.mat'))

crecS = 0; % counter for recording different cases
for ne = 1:Ne
    if (DCS(ne)>0.1)&&(Ncxst(ne)>1)
        tic
%         disp(ne)
        crecS = crecS+1;
        N = NS(ne);
        Nc = Ncxst(ne);
        indx = 1:N;
        SSindx = indx(CompS(1:N,ne)>1e-6);
        pHrS(crecS) = pH(ne);
        
        cip = ciSS(SSindx,SSindx,ne);
        cntSpeciesS(crecS,1:Ns0) = hist(Cxst(1:Ncxst(ne),ne),1:Ns0);
        cntIntTypeS(crecS,1) = sum(sum(cip<0))-Nc; % number of negative interaction coefficients
        cntIntTypeS(crecS,2) = sum(sum(cip>0)); % number of positive interaction coefficients
        CompDevS(crecS) = max(DCd(ne,:));
        CompDevAllS(crecS,1:Nd) = DCd(ne,:);
        pHdCaseS(crecS) = max(pHdd(ne,:));
        pHdCaseAllS(crecS,1:Nd) = pHdd(ne,:);
        
    end
end

FracPosIntS = cntIntTypeS(:,2)./sum(cntIntTypeS,2);
figure
plot(pHdCaseS(FracPosIntS>0),CompDevS(FracPosIntS>0),'.')
xlabel('pH fluctuation amplitude')
ylabel('Composition deviation')

NpHdbin = 10;
pHdbinS = linspace(min(pHdCaseS),max(pHdCaseS),NpHdbin+1);
pHdCaseAllSlin = reshape(pHdCaseAllS,1,crecS*Nd);
CompDevAllSlin = reshape(CompDevAllS,1,crecS*Nd);
GbinS = zeros(1,crecS*Nd);
for cpb = 1:NpHdbin
    selbin = (pHdCaseAllSlin>pHdbinS(cpb))&(pHdCaseAllSlin<=pHdbinS(cpb+1));
    DevOverpHSm(cpb) = mean(CompDevAllSlin(selbin));
    DevOverpHSsd(cpb) = std(CompDevAllSlin(selbin));
    GbinS(selbin) = cpb;
end
figure
errorbar(pHdbinS(1:NpHdbin)+(pHdbinS(2)-pHdbinS(1))/2,DevOverpHSm,DevOverpHSsd)
xlabel('pH fluctuation amplitude')
ylabel('Composition deviation')

figure
boxplot(CompDevAllSlin,GbinS,'positions',[0 pHdbinS(1:NpHdbin)+(pHdbinS(2)-pHdbin(1))/2+0.05],'outliersize',0.1,'width',0.02)
xlabel('pH fluctuation amplitude')
ylabel('Composition deviation')
set(gca,'XTick',0:0.2:1.2,'XTickLabels',0:0.2:1.2)
xlim([0 1.2])

for cb = 1:crecS
    RNDindx = ceil(Ns0*rand(1,Nsp));
    cir = ci(RNDindx,RNDindx);
    cntIntType0S(cb,1) = sum(sum(cir<0))-Nsp; % number of negative interaction coefficients
    cntIntType0S(cb,2) = sum(sum(cir>0)); % number of positive interaction coefficients
end
FracPosInt0S = cntIntType0S(:,2)./sum(cntIntType0S,2);
figure
plot(FracPosInt0S,FracPosIntS,'.')
hold on
plot([0 0.5],[0 0.5],'r')
xlabel('Initial fraction of facilitation')
ylabel('Final fraction of facilitation')
xlim([0 0.5])
ylim([0 0.5])

figure
bar(sum(cntSpeciesS))
xlabel('Species')
ylabel('Persistance Frequency')

NpHbin = 10;
pHbinS = linspace(min(pHrS),max(pHrS),NpHbin+1);
for cpb = 1:NpHbin
    selbin = (pHrS>pHbinS(cpb))&(pHrS<=pHbinS(cpb+1));
    SpeciesFreqS(:,cpb) = sum(cntSpeciesS(selbin,:));
end
figure
semilogy(pHbinS(1:NpHbin)+(pHbinS(2)-pHbinS(1))/2,1/crecS*SpeciesFreqS')
xlabel('pH')
ylabel('Frequency of persistance')
legend('1821','1828','1850','1867','1989','1839')


[hR, cR] = hist(CompDev,logspace(-2,0,8));
[hS, cS] = hist(CompDevS,logspace(-2,0,8));
figure
errorbar(cR,1/sum(hR)*hR,1/sum(hR)*sqrt(hR),'bo-')
hold on
errorbar(cS,1/sum(hS)*hS,1/sum(hS)*sqrt(hS),'rs:')
set(gca,'YScale','log','XScale','log')
xlabel('Dissimilarity with pH fluctuations')
ylabel('Probability distribution')
legend('Robust communities','Sensitive communities')
ylim([8e-3 1])

save(strcat('CompareRobustSensitive_',infile,'.mat'))
