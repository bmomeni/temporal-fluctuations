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

NOL = 0.1; % niche overlap factor

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

rS(3,:) = [1.287025	1.33515	1.47155	1.547075	1.5192	1.4913	1.7351	1.3384	1.20435];
KS(3,:) = [0.4185	0.3125	0.21375	0.3555	0.3295	0.3295	0.34	0.3135	0.234666667];

rS(4,:) = [1.15395	1.321775	1.454425	1.62205	1.56675	1.665875	1.566825	1.373525	1.249525];
KS(4,:) = [0.33125	0.194	0.18325	0.347	0.31325	0.2865	0.22875	0.19975	0.15925];

rS(5,:) = [1.179725,1.341775,1.393475,1.59275,1.5304,1.544475,1.625725,1.3022,1.235875];
KS(5,:) = [0.3175,0.19325,0.18175,0.333,0.3155,0.292,0.28975,0.27875,0.2325];

rS(6,:) = [1.214675	1.4241	1.50545	1.551133333	1.56515	1.61765	1.492225	1.31545	1.119];
KS(6,:) = [0.33175	0.2315	0.201	0.318	0.31375	0.297	0.25525	0.2325	0.2105];

%% Interaction coefficients
ci = -NOL*ones(Ns0,Ns0) + (NOL-1)*eye(Ns0,Ns0);

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
    
    % Niche overlap
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
save(strcat('EnsembleRobustCoexistenceCS_NOL',num2str(NOL),'_Rnd0LVsa_BC_Ngen',num2str(Ngen),'_fp',num2str(fp*100),'_Ne',num2str(Ne),'_rndseed',num2str(rndseed0),'.mat'))
