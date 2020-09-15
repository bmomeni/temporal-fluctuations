%% Finding examples of robust coexistence among nasal microbes
% Ne cases are tested, in each case, a subset of species are put together,
% simulated for 100 generations, and if there is coexistence, a range of
% dilution rates (d/1.5 to d*1.5) are tested to confrim that the
% coexistence is robust and not sensitive to growth conditions.

% "Exp": based on experimental data
% "Spt0LVsa": LV parameters estimated based on supernatant assays;
% populations kept steady (no decline) after reaching stationary phase
% adaptive time-step; diagonal terms of ci set to -1

clear

rndseed0 = 7239;
rng(rndseed0,'twister');

%% Simulation parameters
Nsp = 20; % number of species/strains in each assembly
Ns0 = 6; % number of characteristic strains sampled from
Ne = 10000; % number of cases tested
Ngen = 100; % total number of generations of growth is 100
pHrng = 5.1:0.3:7.5; % Experimental pH range
dmin = 0.03; % min dilution factor
dmax = 0.3; % max dilution factor
dt = 0.02; % simulation time-step
fp = 0.2; % fluctuation parameter when sampling from a type species

%% pH response of species tested
% Species 1 (1821)
% Species 2 (1828)
% Species X (1996) * not used with supernatant
% Species 3 (1850)
% Species 4 (1867)
% Species 5 (1989)
% Species 6 (1839)
rS(1,:) = [0,0.5206,1.4365,1.26905,0.862566667,0.947225,0.991133333,1.039225,0.99015];
KS(1,:) = [0.00075,0.64125,0.6565,0.69,0.65925,0.75225,0.7575,0.8152,0.86];

rS(2,:) = [1.98195,2.0067,2.186775,2.271375,2.149033333,2.145,1.9969,2.047575,2.184925];
KS(2,:) = [0.35675,0.362,0.46175,0.50125,0.5565,0.5765,0.59,0.45275,0.45525];

% rS(4,:) = [0.825025	1.0873	1.146325	1.314125	1.468425	1.497225	1.457475	1.43555	1.464725];
% KS(4,:) = [0.3705	0.42725	0.393	0.351	0.41025	0.417	0.45775	0.441	0.425];

rS(3,:) = [1.287025	1.33515	1.47155	1.547075	1.5192	1.4913	1.7351	1.3384	1.20435];
KS(3,:) = [0.4185	0.3125	0.21375	0.3555	0.3295	0.3295	0.34	0.3135	0.234666667];

rS(4,:) = [1.15395	1.321775	1.454425	1.62205	1.56675	1.665875	1.566825	1.373525	1.249525];
KS(4,:) = [0.33125	0.194	0.18325	0.347	0.31325	0.2865	0.22875	0.19975	0.15925];

rS(5,:) = [1.179725,1.341775,1.393475,1.59275,1.5304,1.544475,1.625725,1.3022,1.235875];
KS(5,:) = [0.3175,0.19325,0.18175,0.333,0.3155,0.292,0.28975,0.27875,0.2325];

rS(6,:) = [1.214675	1.4241	1.50545	1.551133333	1.56515	1.61765	1.492225	1.31545	1.119];
KS(6,:) = [0.33175	0.2315	0.201	0.318	0.31375	0.297	0.25525	0.2325	0.2105];

%% Data from 'Compiled data.xlsx
% Growth in supernatant of other species at ph 7.2
% Species 1 (1821)
% Species 2 (1828)
% Species X (1838) * not studied with pH
% Species Y (2000) * not studied with pH
% Species 3 (1850)
% Species 4 (1867)
% Species 5 (1989)
% Species 6 (1839)
rSPT = [0.9765	0.1099	0.3584	0.1774	0.2312	0	0.9084	0.872	0.52;
    2.574	0.5159	0.2467	1.968	0.6376	1.190	2.979	2.237	0.5978;
    2.031   0.97    1.187   2.785   1.564   2.021   3.489   1.822   1.730;
    0.3920	0.0836	0.07908	0.2994	0.1930	0.3168	0.2946	0.439	0.1602;
    1.428	0.08458	0.1167	0.4812	0.01593	0.3351	1.653	1.079	0;
    1.509	0.05878	0.00458	0.5562	0       0.3279	1.740	1.006	0;
    0.2918	0       0.02208	0.2316	0.08858	0.1387	0.1720	0.2022	0.0796;
    1.454	0.08708	0.249	0.5034	0.04922	0.1827	1.725	0.9412	0];

KSPT = [0.7089	0.2022	0.02792	0.1282	0.189	0.00342	0.5739	0.712	0.4022;
    0.6527	0.1747	0.1829	0.4522	0.1499	0.1772	0.4704	0.4294	0.1552;
    0.642	0.2522	0.2787	0.4012	0.2442	0.2707	0.3519	0.5377	0.2342;
    0.2654	0.01467	0.01117	0.2222	0.00342	0.02942	0.2739	0.1562	0.000917;
    0.2072	0.01892	0.01767	0.1999	0.03117	0.03067	0.2539	0.1112	0.001167;
    0.3422	0.1009	0.1289	0.1447	0.1639	0.1704	0.3482	0.2154	0.1734;
    0.2862	0.00217	0.01442	0.1804	0.05192	0.07467	0.3062	1.1022	0.05392;
    0.2494	0.04467	0.02117	0.1519	0.01217	0.03317	0.2402	0.1947	0.00642];

seli = [1 2 5 6 7 8]; % selected strains from the panel
Nsel = length(seli); % # of selected strains from the panel
r0 = rSPT(seli,1);
K0 = KSPT(seli,1);

% Rearranging data (first column is fresh media)
rSPTn = rSPT(seli,seli+1); % growth rates in supernatants
KSPTn = KSPT(seli,seli+1); % carrying capacity in supernatants

%% Interaction coefficients
ci = (KSPTn - K0*ones(1,Nsel))./(ones(Nsel,1)*K0');
% set diagonal terms to -1
ci = ci - diag(diag(ci)) - eye(Ns0,Ns0);

S0i = 1e-4; % average initial density of each population (in OD)

%% Record of cases tested
d0 = zeros(1,Ne); % dilution factor
pH = zeros(1,Ne); % pH
NS = zeros(1,Ne); % number of species picked initially
Ncxst = zeros(1,Ne); % number of species that reached coexistence
DCS = zeros(1,Ne); % deviation in composition from the reference
NsmplS = zeros(Nsp,Ne); % record of type species randomly modified and sampled
CompS = zeros(Nsp,Ne); % final composition after coexistence
rSS = zeros(Nsp,Ne); % growth rates of species
KSS = zeros(Nsp,Ne); % carrying capacities of species
ciSS = zeros(Nsp,Nsp,Ne); % record of instances of interaction coefficients

tic
%% Surveying Ne cases
ND = 3;
indx = 1:Nsp;
for ne = 1:Ne
    if mod(ne,50)==1
        disp(ne)
        toc
        tic
    end
    %     tic
    d0(ne) = exp(log(dmin)+rand(1)*(log(dmax)-log(dmin)));
    drng = d0(ne)*[1/1.5, 1, 1.5]; % dilution rate at each transfer
    
    pH(ne) = min(pHrng)+rand(1)*(max(pHrng)-min(pHrng)); % test pH
    
    N = 1+ceil(rand(1)*(Nsp-1)); % number of species kept in this instance
    
    r = zeros(N,1);
    K = zeros(N,1);
    
    % Which type species each sampled strain belongs to
    Nsmpl = 1+floor(rand(N,1)*Ns0);
    Nsmpl(Nsmpl>Ns0) = Ns0;
    
    % interaction coefficients
    cip = ci(Nsmpl,Nsmpl).*((1-fp)+2*fp*rand(N,N));
    
    % initial population density (cells/ml)
    S0 = S0i*rand(N,1);
    
    for n = 1:N
        % parameters at given pH
        r(n) = interp1(pHrng,rS(Nsmpl(n),:),pH(ne))*((1-fp)+2*fp*rand(1)); % basal growth rates of species 1
        K(n) = interp1(pHrng,KS(Nsmpl(n),:),pH(ne))*((1-fp)+2*fp*rand(1)); % basal yields of species 1
    end
    d = d0(ne);
    
    S = S0;
    Sd = zeros(size(S));
    ct = 0;
    td = 0;
    Gen = 0;
    while (Gen < Ngen) % duration of simulations
        ct = ct+1;
        
        % assuming logistic growth
        re = (r.*(1+1./K.*(cip*S)).*((1+1./K.*(cip*S))>0)).*(S>0);
        Su = S + dt*(re-d).*S; % updated pop. sizes because of growth
        td = td+dt;
        
        Gen = Gen + d*dt/log(2);
        dt = 0.1/max(abs(re-d));
        S = Su;
    end
    
    Comp0 = 1/sum(S)*S;
    Comp0(Comp0<1e-6) = 0;
    Ncxst(ne) = sum(Comp0 > 0);
    
    %% Test robustness against dilution rate
    DC = zeros(ND,1);
    if Ncxst(ne)>1
        cD = 0;
        for d = drng
            cD = cD+1;
            
            S = S0i*Comp0;
            
            td = 0;
            cntd = 0;
            cnt = 0;
            Gen = 0;
            while (Gen < Ngen) % duration of simulations
                cnt = cnt+1;
                
                % assuming logistic growth within each well
                re = (r.*(1+1./K.*(cip*S)).*((1+1./K.*(cip*S))>0)).*(S>0);
                Su = S + dt*(re-d).*S; % updated pop. sizes because of growth
                td = td+dt;
                
                Gen = Gen + d*dt/log(2);
                dt = 0.1/max(abs(re-d));
                S = Su;
            end
            Compf = 1/sum(S)*S;
            % Bray-Curtis measure of composition dissimmilarity
            CompDist = f_dis([Compf,Comp0]','BC');
            DC(cD) = CompDist(1,2);
        end
        DCS(ne) = max(DC);
    end
    
    NS(ne) = N; % number of species at inoculation
    NsmplS(1:N,ne) = Nsmpl; % type strains sampled in this instance
    CompS(1:N,ne) = Comp0(1:N); % final composition of assembled community
    rSS(1:N,ne) = r; % growth rate of strains at given pH
    KSS(1:N,ne) = K; % carrying capacities of strains at given pH
    ciSS(1:N,1:N,ne) = cip; % interaction coefficients
end


%% Recording the results
save(strcat('EnsembleRobustCoexistenceCS_RndExpSpt0LVsa_BC_Ngen',num2str(Ngen),'_fp',num2str(fp*100),'_Ne',num2str(Ne),'_rndseed',num2str(rndseed0),'.mat'))
