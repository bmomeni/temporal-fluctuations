clear

load('pHCheckC_fpH0.2_EnsembleRobustCoexistenceCS_RndExpSpt0LVsa_BC_Ngen100_fp20_Ne10000_rndseed7239.mat');

pHdd = zeros(Ne,Nd);
cim = 0;
for ne = 1:Ne
    if (DCS(ne)<0.1)&&(Ncxst(ne)>1)&&(DCd(ne,1)<0.1)
        cim = cim+1;
        pH0 = pH(ne);
        pHdm = min(pH(ne)-min(pHrng),max(pHrng)-pH(ne));
        pHdd(ne,:) = pHdm*linspace(0,1,Nd);
        pHim = pHdm;
        DCim(cim,:) = DCd(ne,:);
    else
        pHdd(ne,:) = NaN;
    end
end

% figure
% plot(pHdd,DCd,'k.')
% xlabel('pH fluctuation')
% ylabel('Composition deviation')
% xlim([0 max(max(pHdd))])
% ylim([0 1])

Nn = Nd*sum((DCS<0.1).*(Ncxst>1));
disp('Number of cases examined:')
disp(Nn)

Nch = sum(sum(DCd>0.1));
disp('Number of significant composition deviations:')
disp(Nch)

disp('Freaction of cases affected by pH fluctuations')
disp(Nch/Nn)

pHddl = reshape(pHdd,1,Ne*Nd);
DCdl = reshape(DCd,1,Ne*Nd);
cps = 0;
pHsrng = 0.1:0.1:1;
for pHstep = pHsrng
    cps = cps+1;
    Nn(cps) = sum((pHddl<pHstep)&(pHddl>=pHstep-0.1));
    Nch2(cps) = sum(DCdl((pHddl<pHstep)&(pHddl>=pHstep-0.1))>0.2);
    Nch1(cps) = sum(DCdl((pHddl<pHstep)&(pHddl>=pHstep-0.1))>0.1) - Nch2(cps);

    Nch0(cps) = Nn(cps) - Nch1(cps) - Nch2(cps);

    [phi,pci] = binofit(cumsum([Nch0(cps); Nch1(cps); Nch2(cps)])',Nn(cps)*[1 1 1]);
    FchLCI(:,cps) = pci(:,1);
    FchHCI(:,cps) = pci(:,2); 
end

Fch = 1./(ones(3,1)*Nn).*[Nch0; Nch1; Nch2];
FchC = cumsum(Fch);

figure
bar(pHsrng-0.05,Fch',1,'stacked')
hold on
for ct = 1:3
    errorbar(pHsrng-0.05,FchC(ct,:),FchLCI(ct,:)-FchC(ct,:),FchHCI(ct,:)-FchC(ct,:),'k.')
end
colormap([0.2 0 1; 0.7 0.7 0.9; 1 0.6 0])
xlabel('pH fluctuation amplitude')
ylabel('Occurrence')
legend('<0.1','>0.1, <0.2','>0.2')
ylim([0 1])
% 
% figure
% imagesc(DCim)