clear

mcirng = 0.1:0.1:0.8;
cmci = 0;
Nmci = length(mcirng);

pHdt = 0.5; % targeted pH fluctuation
pHdtd = 0.1; % acceptable range around the targeted value

for mcid = mcirng
    cmci = cmci+1;
    
    load(strcat('pHCheck_EnsembleRobustCoexistenceCS_mPDI',num2str(10*mcid),'_RndExpSpt0LVsa_BC_Ngen100_fp20_Ne10000_rndseed7239.mat'));
    
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
    
%     Nn = Nd*sum((DCS<0.1).*(Ncxst>1));
%     disp('Number of cases examined:')
%     disp(Nn)
%     
%     Nch = sum(sum(DCd>0.1));
%     disp('Number of composition deviations >0.1:')
%     disp(Nch)
%     
%     disp('Freaction of cases affected by pH fluctuations')
%     disp(Nch/Nn)
    
    pHddl = reshape(pHdd,1,Ne*Nd);
    DCdl = reshape(DCd,1,Ne*Nd);
    
    pHInRange = abs(pHddl-pHdt)<pHdtd;
    Nn(cmci) = sum(pHInRange);
    Nch2(cmci) = sum(DCdl(pHInRange)>0.2);
    Nch1(cmci) = sum(DCdl(pHInRange)>0.1) - Nch2(cmci);
    
    Nch0(cmci) = Nn(cmci) - Nch1(cmci) - Nch2(cmci);
    
    [phi,pci] = binofit(cumsum([Nch0(cmci); Nch1(cmci); Nch2(cmci)])',Nn(cmci)*[1 1 1]);
    FchLCI(:,cmci) = pci(:,1);
    FchHCI(:,cmci) = pci(:,2);
end

Fch = 1./(ones(3,1)*Nn).*[Nch0; Nch1; Nch2];
FchC = cumsum(Fch);

figure
hold on
bar(mcirng,Fch',0.2,'stacked')
for ct = 1:3
    errorbar(mcirng,FchC(ct,:),FchLCI(ct,:)-FchC(ct,:),FchHCI(ct,:)-FchC(ct,:),'k.')
end
xlabel('Degree of interaction change by pH')
ylabel('Occurrence')
legend('<0.1','>0.1, <0.2','>0.2')
colormap([0.2 0 1; 0.7 0.7 0.9; 1 0.6 0])
ylim([0 1])
xlim([0 0.9])
