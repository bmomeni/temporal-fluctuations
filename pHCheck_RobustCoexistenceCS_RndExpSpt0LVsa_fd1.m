clear

infile = 'EnsembleRobustCoexistenceCS_RndExpSpt0LVsa_BC_Ngen100_fp20_Ne10000_rndseed7239';
load(strcat(infile,'.mat'))


Nd = 8; % number of pH frequencies tested
fpHrng = logspace(-1.3,0.3,Nd);
pHd = 0.1; % frequency of pH variation (1/hr); higher = faster
dt = 0.1;

Compfd = zeros(Nsp,Ne,Nd);
DCd = zeros(Ne,Nd);
fpHd = zeros(Ne,Nd);
Cxst = zeros(Nsp,Ne);

for ne = 1:Ne
    if (DCS(ne)<0.1)&&(Ncxst(ne)>1)&&(abs(pH(ne)-mean(pHrng))<(max(pHrng)-mean(pHrng)-pHd))
        tic
        disp(ne)
        N = NS(ne);
        Nc = Ncxst(ne);
        indx = 1:N;
        SSindx = indx(CompS(1:N,ne)>1e-6);
        pH0 = pH(ne);
        d = d0(ne);
        
        r = rSS(SSindx,ne);
        K = KSS(SSindx,ne);
        Cxst(1:Ncxst(ne),ne) = NsmplS(SSindx,ne);
        cip = ciSS(SSindx,SSindx,ne);
        
        nd = 0;
        for fpH = fpHrng
            nd = nd+1;
            td = 0;
            tu = round(0.05/fpH);
            tph = 0;
            Gen = 0;
            cnt = 0;
            % initial population density (cells/ml)
            S0 = S0i*CompS(SSindx,ne);
            S = S0;
            while (Gen < Ngen) % between dilution steps
                cnt = cnt+1;
                
                tph = tph+dt;
                % find parameters for current pH
                if tph > tu
                    tph = 0;
                    pHt = pH0+pHd*sin(2*pi*fpH*td);% parameters at given pH
                    for n = 1:Ncxst(ne)
                        % parameters at given pH
                        r(n) = rSS(SSindx(n),ne)/interp1(pHrng,rS(Cxst(n,ne),:),pH0)*interp1(pHrng,rS(Cxst(n,ne),:),pHt); % basal growth rates of species 1
                        K(n) = KSS(SSindx(n),ne)/interp1(pHrng,KS(Cxst(n,ne),:),pH0)*interp1(pHrng,KS(Cxst(n,ne),:),pHt); % basal yields of species 1
                        if K(n)<1e-4
                            K(n) = 1e-4;
                        end
                    end
                end
                % assuming logistic growth within each well
                re = (r.*(1+1./K.*(cip*S)).*((1+1./K.*(cip*S))>0)).*(S>0);
                Su = S + dt*(re-d).*S; % updated pop. sizes because of growth
                td = td+dt;
                
                Gen = Gen + d*dt/log(2);
                dt = min(0.2,0.1/max(abs(re-d)));
                S = Su;
                
            end
            tcrng = linspace(td,td+1/fpH-mod(td,1/fpH),ceil(3*(1/fpH-mod(td,1/fpH))/dt));
            dtc = (1/fpH-td)/(ceil(3*(1/fpH-td)/dt)-1);
            for td = tcrng % complete the cycle for a round number periods
                cnt = cnt+1;
                
                tph = tph+dtc;
                % find parameters for current pH
                if tph > tu
                    tph = 0;
                    pHt = pH0+pHd*sin(2*pi*fpH*td);% parameters at given pH
                    for n = 1:Ncxst(ne)
                        % parameters at given pH
                        r(n) = rSS(SSindx(n),ne)/interp1(pHrng,rS(Cxst(n,ne),:),pH0)*interp1(pHrng,rS(Cxst(n,ne),:),pHt); % basal growth rates of species 1
                        K(n) = KSS(SSindx(n),ne)/interp1(pHrng,KS(Cxst(n,ne),:),pH0)*interp1(pHrng,KS(Cxst(n,ne),:),pHt); % basal yields of species 1
                        if K(n)<1e-4
                            K(n) = 1e-4;
                        end
                    end
                end
                % assuming logistic growth within each well
                re = (r.*(1+1./K.*(cip*S)).*((1+1./K.*(cip*S))>0)).*(S>0);
                Su = S + dt*(re-d).*S; % updated pop. sizes because of growth
                
                Gen = Gen + d*dt/log(2);
                dt = min(0.2,0.1/max(abs(re-d)));
                S = Su;
                
            end
            Compfd(1:Ncxst(ne),ne,nd) = 1/sum(S)*S;
            CompRef0 = CompS(SSindx,ne);
            CompRef = CompRef0(CompRef0>1e-6);
            % Bray-Curtis measure of composition dissimmilarity
            CompDist = f_dis([Compfd(1:Ncxst(ne),ne,nd),CompRef]','BC');
            DCd(ne,nd) = CompDist(1,2);
            fpHd(ne,nd) = pHd;
            
        end
        disp(DCd(ne,:))
        toc
    end
end

save(strcat('pHCheck_fpHrng_pHd',num2str(pHd),'_',infile,'.mat'))
