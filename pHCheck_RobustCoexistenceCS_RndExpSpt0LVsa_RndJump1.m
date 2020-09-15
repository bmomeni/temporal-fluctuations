clear

%% Need to update: preset r and K at p+/-delta_p
% assign time between transitions
% switch between parameters when Tt has lapsed
%

infile = 'EnsembleRobustCoexistenceCS_RndExpSpt0LVsa_BC_Ngen100_fp20_Ne10000_rndseed7239';
load(strcat(infile,'.mat'))


Nd = 11; % number of pH frequencies tested
ftrng = [0.0187    0.0260    0.0361 logspace(-1.3,0.3,Nd-3)];
pHd = 0.5; % amplitude of pH variation
dt = 0.1;

Compfd = zeros(Nsp,Ne,Nd);
DCd = zeros(Ne,Nd);
fpHd = zeros(Ne,Nd);
Cxst = zeros(Nsp,Ne);

for ne = 1:Ne
    if (DCS(ne)<0.1)&&(Ncxst(ne)>1)&&(abs(pH(ne)-mean(pHrng))<(max(pHrng)-mean(pHrng)-pHd))
        disp(ne)
        N = NS(ne);
        Nc = Ncxst(ne);
        indx = 1:N;
        SSindx = indx(CompS(1:N,ne)>1e-6);
        pH0 = pH(ne);
        d = d0(ne);
        Cxst(1:Ncxst(ne),ne) = NsmplS(SSindx,ne);
        rt = zeros(Nc,2);
        Kt = zeros(Nc,2);
        for n = 1:Ncxst(ne)
            rt(n,1) = rSS(SSindx(n),ne)/interp1(pHrng,rS(Cxst(n,ne),:),pH0)*interp1(pHrng,rS(Cxst(n,ne),:),pH0-pHd); % basal growth rates at pH1
            Kt(n,1) = KSS(SSindx(n),ne)/interp1(pHrng,KS(Cxst(n,ne),:),pH0)*interp1(pHrng,KS(Cxst(n,ne),:),pH0-pHd); % basal yields at pH1
            if Kt(n,1)<1e-4
                Kt(n,1) = 1e-4;
            end
            rt(n,2) = rSS(SSindx(n),ne)/interp1(pHrng,rS(Cxst(n,ne),:),pH0)*interp1(pHrng,rS(Cxst(n,ne),:),pH0+pHd); % basal growth rates at pH2
            Kt(n,2) = KSS(SSindx(n),ne)/interp1(pHrng,KS(Cxst(n,ne),:),pH0)*interp1(pHrng,KS(Cxst(n,ne),:),pH0+pHd); % basal yields at pH2
            if Kt(n,2)<1e-4
                Kt(n,2) = 1e-4;
            end
            rt(n,3) = rSS(SSindx(n),ne); % basal growth rates at pH0
            Kt(n,3) = KSS(SSindx(n),ne); % basal yields at pH0
            if Kt(n,3)<1e-4
                Kt(n,3) = 1e-4;
            end
        end
        cip = ciSS(SSindx,SSindx,ne);
        
        nd = 0;
        for ft = ftrng
            nd = nd+1;
            td = 0;
            tph = 0;
            Gen = 0;
            cnt = 0;
            ptt = 1;
            % initial population density (cells/ml)
            S0 = S0i*CompS(SSindx,ne);
            S = S0;
            while (Gen < Ngen) % between dilution steps
                
                Tt = -1/ft*log(1-rand(1)); % time until next transition
                ttrng = linspace(0,Tt,ceil(Tt/dt)+1);
                dt = ttrng(2)-ttrng(1);
                for tt = ttrng
                    cnt = cnt+1;
                    
                    % assuming logistic growth within each well
                    re = (rt(:,ptt).*(1+1./Kt(:,ptt).*(cip*S)).*((1+1./Kt(:,ptt).*(cip*S))>0)).*(S>0);
                    Su = S + dt*(re-d).*S; % updated pop. sizes because of growth
                    td = td+dt;
                    
                    Gen = Gen + d*dt/log(2);
                    ptt = 3 - ptt; % switch from 2 to 1 or from 1 to 2
                end
                dt = min(0.2,0.1/max(abs(re-d)));
                S = Su;
                
            end
            % add a pH0 cycle to finish all runs on the same cycle
            Tt = 1/ft; % time until next transition
            ttrng = linspace(0,Tt,ceil(Tt/dt)+1);
            dt = ttrng(2)-ttrng(1);
            for tt = ttrng
                cnt = cnt+1;
                
                % assuming logistic growth within each well
                re = (rt(:,3).*(1+1./Kt(:,3).*(cip*S)).*((1+1./Kt(:,3).*(cip*S))>0)).*(S>0);
                Su = S + dt*(re-d).*S; % updated pop. sizes because of growth
                td = td+dt;
                
                Gen = Gen + d*dt/log(2);
            end
            S = Su;
            
            Compfd(1:Ncxst(ne),ne,nd) = 1/sum(S)*S;
            CompRef0 = CompS(SSindx,ne);
            CompRef = CompRef0(CompRef0>1e-6);
            % Bray-Curtis measure of composition dissimmilarity
            CompDist = f_dis([Compfd(1:Ncxst(ne),ne,nd),CompRef]','BC');
            DCd(ne,nd) = CompDist(1,2);
            fpHd(ne,nd) = pHd;
            
        end
    end
end

save(strcat('pHCheck_RndJump_ftrng_pHd',num2str(pHd),'_',infile,'.mat'))
