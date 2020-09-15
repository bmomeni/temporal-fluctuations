clear

load('pHCheck_RndJump_ftrng_pHd0.25_EnsembleRobustCoexistenceCS_RndExpSpt0LVsa_BC_Ngen100_fp20_Ne10000_rndseed7239.mat');
% load('pHCheck_RndJump_ftrng_pHd0.5_EnsembleRobustCoexistenceCS_RndExpSpt0LVsa_BC_Ngen100_fp20_Ne10000_rndseed7239.mat');

% figure
% plot((ones(Ne,1)*fpHrng).*(0.9+0.2*rand(Ne,8)),DCd,'k.')
% xlabel('pH fluctuation freq. (1/hr)')
% ylabel('Composition deviation')
% set(gca,'XScale','log')

cfp = 0;
Nft = length(ftrng);

for ft = ftrng
    cfp = cfp+1;
    selindx = (DCS<0.1)&(Ncxst>1)&(abs(pH-mean(pHrng))<(max(pHrng)-mean(pHrng)-pHd));
    Nn(cfp) = sum(selindx);
    Nch2(cfp) = sum(DCd(selindx,cfp)>0.2);
    Nch1(cfp) = sum(DCd(selindx,cfp)>0.1) - Nch2(cfp);
    Nch0(cfp) = Nn(cfp) - Nch1(cfp) - Nch2(cfp);

    [phi,pci] = binofit(cumsum([Nch0(cfp); Nch1(cfp); Nch2(cfp)])',Nn(cfp)*[1 1 1]);
    FchLCI(:,cfp) = pci(:,1);
    FchHCI(:,cfp) = pci(:,2); 
end

Fch = 1./(ones(3,1)*Nn).*[Nch0; Nch1; Nch2];
FchC = cumsum(Fch);
figure
hold on
for cfp = 1:Nft
    bar([0 ftrng(cfp)],[zeros(3,1) Fch(:,cfp)]',0.1*ftrng(cfp)^0.03,'stacked')
end
for ct = 1:3
    errorbar(ftrng,FchC(ct,:),FchLCI(ct,:)-FchC(ct,:),FchHCI(ct,:)-FchC(ct,:),'k.')
end
xlabel('pH transition freq. (1/hr)')
ylabel('Occurrence')
legend('<0.1','>0.1, <0.2','>0.2')
colormap([0.2 0 1; 0.7 0.7 0.9; 1 0.6 0])
ylim([0 1])
xlim([0.01 3])
set(gca,'XScale','log')

