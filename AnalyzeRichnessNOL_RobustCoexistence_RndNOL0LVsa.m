clear

NOLrng = 0.1:0.1:1;
nNOL = 0;
for NOLd = NOLrng
    nNOL = nNOL+1;
    load(strcat('EnsembleRobustCoexistenceCS_NOL',num2str(NOLd),'_Rnd0LVsa_BC_Ngen100_fp20_Ne10000_rndseed7239.mat'))
    
    Cxst = zeros(Nsp,Ne);
    crec = 0; % counter for recording different cases
    pHd = 0.5;
    for ne = 1:Ne
        if (DCS(ne)<0.1)&&(Ncxst(ne)>1)&&(abs(pH(ne)-mean(pHrng))<(max(pHrng)-mean(pHrng)-pHd))
            crec = crec+1;
            N = NS(ne);
            Nc = Ncxst(ne);
            indx = 1:N;
            SSindx = indx(CompS(1:N,ne)>1e-6);
            InitialSpecies = NsmplS(1:N,ne);
            FinalSpecies = InitialSpecies(CompS(1:N,ne)>1e-6);
            Richness(crec) = length(unique(FinalSpecies));
            Cxst(1:Nc,ne) = NsmplS(SSindx,ne);
            pHr(crec) = pH(ne);
            
            cir = ci(InitialSpecies,InitialSpecies);
            cntIntType0(crec,1) = sum(sum(cir<0))-N; % number of negative interaction coefficients
            cntIntType0(crec,2) = sum(sum(cir>0)); % number of positive interaction coefficients
            
            cip = ci(FinalSpecies,FinalSpecies);
            cntSpecies(crec,1:Ns0) = hist(Cxst(1:Ncxst(ne),ne),1:Ns0);
            cntIntType(crec,1) = sum(sum(cip<0))-Nc; % number of negative interaction coefficients
            cntIntType(crec,2) = sum(sum(cip>0)); % number of positive interaction coefficients
        end
    end
    Nn(nNOL) = crec;
    MeanRichness(nNOL) = mean(Richness(1:crec));
    StdRichness(nNOL) = std(Richness(1:crec));
end

figure
errorbar(NOLrng,1/Ne*Nn,1/Ne*sqrt(Nn))
xlabel('Interspecies niche overlap')
ylabel('Likelihood of coexistence')
ylim([0 1])
xlim([0 1.1])

figure
errorbar(NOLrng,MeanRichness,StdRichness)
xlabel('Interspecies niche overlap')
ylabel('Mean species richness')
xlim([0 1.1])
