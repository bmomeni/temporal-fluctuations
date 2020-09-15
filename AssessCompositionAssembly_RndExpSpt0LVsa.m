clear

infile = {
    'pHCheckAssemblyCS_RndExpSpt0LVsa_pHd_Ngen100_fp20_fpH0.02_Ne1000_rndseed7239'
    'pHCheckAssemblyCS_RndExpSpt0LVsa_pHd_Ngen100_fp20_fpH0.04_Ne1000_rndseed7239'
    'pHCheckAssemblyCS_RndExpSpt0LVsa_pHd_Ngen100_fp20_fpH0.08_Ne1000_rndseed7239'
    'pHCheckAssemblyCS_RndExpSpt0LVsa_pHd_Ngen100_fp20_fpH0.16_Ne1000_rndseed7239'
    'pHCheckAssemblyCS_RndExpSpt0LVsa_pHd_Ngen100_fp20_fpH0.32_Ne1000_rndseed7239'
    'pHCheckAssemblyCS_RndExpSpt0LVsa_pHd_Ngen100_fp20_fpH0.64_Ne1000_rndseed7239'
    'pHCheckAssemblyCS_RndExpSpt0LVsa_pHd_Ngen100_fp20_fpH1.28_Ne1000_rndseed7239'
    };

NfpH = length(infile);
fpHrng = [0.02 0.04 0.08 0.16 0.32 0.64 1.28];

NpHbin = 10;
Ne = 1000;
Nd = 10;
Ns0 = 6;
CNC = zeros(Ne,Nd,NfpH);
CollectedSpDrop = zeros(NfpH,NpHbin,Ns0);
CollectedSpAdd = zeros(NfpH,NpHbin,Ns0);
for nfph = 1:NfpH
    load(strcat(infile{nfph},'.mat'))
    crec = 0; % counter for recording different cases
    CNC(:,:,nfph) = -Ncxstd + Ncxst'*ones(1,Nd);
    pHdC(:,:,nfph) = pHdd;
    for ne = 1:Ne
        NSpRichness(ne,nfph) = length(unique(NsmplS(CompS(1:NS(ne),ne)>0,ne)));
        for nd = 1:Nd
            NSpRichnessd(ne,nd,nfph) = length(unique(NsmplS(Compfd(1:NS(ne),ne,nd)>0,ne)));
        end
    end
    CNClin = reshape(-Ncxstd + Ncxst'*ones(1,Nd),1,Ne*Nd);
    pHdClin = reshape(pHdd,1,Ne*Nd);
    CNSplin = reshape(-NSpRichnessd(:,:,nfph) + NSpRichness(:,nfph)*ones(1,Nd),1,Ne*Nd);
    pHbin = linspace(min(pHdClin),max(pHdClin),NpHbin+1);
    for cpb = 1:NpHbin
        selbin = (pHdClin>pHbin(cpb))&(pHdClin<=pHbin(cpb+1));
        cntAll(nfph,cpb) = sum(selbin);
        IncRichness(nfph,cpb) = 1/sum(selbin)*sum(CNClin((CNClin>0)&selbin));
        DecRichness(nfph,cpb) = 1/sum(selbin)*sum(CNClin((CNClin<0)&selbin));
        ChangeRichness(nfph,cpb) = sum(CNClin(selbin));
        ChangeSpRichness(nfph,cpb) = sum(CNSplin(selbin));
        cntIncRichness(nfph,cpb) = sum(CNClin(selbin)>0);
        cntDecRichness(nfph,cpb) = sum(CNClin(selbin)<0);

        for tst = 1:Ne*Nd
            if selbin(tst)==1
                tste = ceil(tst/Nd);
                tstd = tst - (tste-1)*Nd;
                Sp0 = unique(NsmplS(CompS(1:NS(tste),tste)>0,tste));
                Spd = unique(NsmplS(Compfd(1:NS(tste),tste,tstd)>0,tste));
                for sp = 1:Ns0
                    CollectedSpAdd(nfph,cpb,sp) = CollectedSpAdd(nfph,cpb,sp) + ismember(sp,setdiff(Spd,Sp0));
                    CollectedSpDrop(nfph,cpb,sp) = CollectedSpDrop(nfph,cpb,sp) + ismember(sp,setdiff(Sp0,Spd));
                end
            end
        end
    end
end

figure
bar(pHbin(1:NpHbin)+(pHbin(2)-pHbin(1))/2,IncRichness')
hold on
bar(pHbin(1:NpHbin)+(pHbin(2)-pHbin(1))/2,DecRichness')
xlabel('pH fluctuation amplitude')
ylabel('Change in richness')

figure
bar(pHbin(1:NpHbin)+(pHbin(2)-pHbin(1))/2,1/Ne*cntIncRichness')
hold on
bar(pHbin(1:NpHbin)+(pHbin(2)-pHbin(1))/2,-1/Ne*cntDecRichness')
xlabel('pH fluctuation amplitude')
ylabel('Occurrence')
ylim([-0.3 0.3])

figure
bar(pHbin(1:NpHbin)+(pHbin(2)-pHbin(1))/2,(1/Nd/Ne*ChangeRichness)')
xlabel('pH fluctuation amplitude')
ylabel('Average change in richness')

figure
bar(pHbin(1:NpHbin)+(pHbin(2)-pHbin(1))/2,(1/Nd/Ne*ChangeSpRichness)')
xlabel('pH fluctuation amplitude')
ylabel('Average change in species richness')

figure
errorbar(pHbin(1:NpHbin)+(pHbin(2)-pHbin(1))/2,mean(ChangeRichness),std(ChangeRichness),'o-')
hold on
errorbar(pHbin(1:NpHbin)+(pHbin(2)-pHbin(1))/2,mean(ChangeSpRichness),std(ChangeSpRichness),'o-')
xlabel('pH fluctuation amplitude')
ylabel('Change in richness')
xlim([0 1.1])

figure
bar(pHbin(1:NpHbin)+(pHbin(2)-pHbin(1))/2,shiftdim(mean(CollectedSpAdd,1),1)./(mean(cntAll)'*ones(1,Ns0)))
hold on
bar(pHbin(1:NpHbin)+(pHbin(2)-pHbin(1))/2,-shiftdim(mean(CollectedSpDrop,1),1)./(mean(cntAll)'*ones(1,Ns0)))
xlabel('pH fluctuation amplitude')
ylabel('Chance of richness change')
ylim([-0.07 0.07])

figure
bar(pHbin(1:NpHbin)+(pHbin(2)-pHbin(1))/2,sum(shiftdim(mean(CollectedSpAdd-CollectedSpDrop,1),1)./(mean(cntAll)'*ones(1,Ns0)),2))
xlabel('pH fluctuation amplitude')
ylabel('Chance of increased richness')

figure
bar(1:Ns0,1/(Ne*Nd)*sum(shiftdim(mean(CollectedSpDrop,1),1)))
hold on
errorbar(1:Ns0,1/(Ne*Nd)*sum(shiftdim(mean(CollectedSpDrop,1),1)),1/(Ne*Nd)*sqrt(sum(shiftdim(mean(CollectedSpDrop,1),1))),'o')
% mm = mean(1/(Ne*Nd)*sum(shiftdim(mean(CollectedSpDrop,1),1)));
% ms = std(1/(Ne*Nd)*sum(shiftdim(mean(CollectedSpDrop,1),1)));
% plot([0 7],[mm mm],'k')
% plot([0 7],[mm mm]+ms,'k:')
% plot([0 7],[mm mm]-ms,'k:')
xlabel('Species')
ylabel('Chance of extinction')
xlim([0 7])
ylim([0 0.05])

figure
bar(1:Ns0,1/(Ne*Nd)*sum(shiftdim(mean(CollectedSpAdd,1),1)))
hold on
errorbar(1:Ns0,1/(Ne*Nd)*sum(shiftdim(mean(CollectedSpAdd,1),1)),1/(Ne*Nd)*sqrt(sum(shiftdim(mean(CollectedSpAdd,1),1))),'o')
% mm = mean(1/(Ne*Nd)*sum(shiftdim(mean(CollectedSpAdd,1),1)));
% ms = std(1/(Ne*Nd)*sum(shiftdim(mean(CollectedSpAdd,1),1)));
% plot([0 7],[mm mm],'k')
% plot([0 7],[mm mm]+ms,'k:')
% plot([0 7],[mm mm]-ms,'k:')
xlabel('Species')
ylabel('Chance of augmentation')
xlim([0 7])
ylim([0 0.05])

save(strcat('AssessRichnessAssembly_',infile{1},'.mat'))
