clear

NOLrng = 0.1:0.1:1;
cNOL = 0;
NNOL = length(NOLrng);

for NOLd = NOLrng
    cNOL = cNOL+1;
    load(strcat('pHCheck_pHd0.5_fpH0.2_EnsembleRobustCoexistenceCS_NOL',num2str(NOLd),'_Rnd0LVsa_BC_Ngen100_fp20_Ne10000_rndseed7239.mat'));
    selindx = (DCS<0.1)&(Ncxst>1)&(abs(pH-mean(pHrng))<(max(pHrng)-mean(pHrng)-pHd));
    Nn(cNOL) = sum(selindx);
    Nch2(cNOL) = sum(DCd(selindx)>0.2);
    Nch1(cNOL) = sum(DCd(selindx)>0.1) - Nch2(cNOL);
    Nch0(cNOL) = Nn(cNOL) - Nch1(cNOL) - Nch2(cNOL);

    [phi,pci] = binofit(cumsum([Nch0(cNOL); Nch1(cNOL); Nch2(cNOL)])',Nn(cNOL)*[1 1 1]);
    FchLCI(:,cNOL) = pci(:,1);
    FchHCI(:,cNOL) = pci(:,2); 
end

Fch = 1./(ones(3,1)*Nn).*[Nch0; Nch1; Nch2];
FchC = cumsum(Fch);
figure
hold on
bar(NOLrng,Fch',0.2,'stacked')
for ct = 1:3
    errorbar(NOLrng,FchC(ct,:),FchLCI(ct,:)-FchC(ct,:),FchHCI(ct,:)-FchC(ct,:),'k.')
end
xlabel('Interspecies niche overlap')
ylabel('Occurrence')
legend('<0.1','>0.1, <0.2','>0.2')
colormap([0.2 0 1; 0.7 0.7 0.9; 1 0.6 0])
ylim([0 1])
xlim([0 1.1])

