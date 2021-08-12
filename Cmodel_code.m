% clear all
% close all
lcolor = 'r';
SUP = 0;
if SUP == 1
    load SpUpvalues2
end
load param
% param = param2;
% load SpUp
%% Model profile depth parameters
D = 0.5 ; % use for when only looking at top 1 meter 
dz = 0.02; %The size of each cell in m (0.01 m = 1 cm)
z = (0:dz:D)';  % Depth vector 
Nz = D/dz+1; %Number of nodes in the z direction (depth)
Aabs = 10e-12; 
An = 6.0221409e23; % Avogadros number
Rpdb = 0.0112372; % In per mil, PDB standard ratio of 13C/12C

% endtime = 10000; % Duration of run in years
endtime = 1000;
end_date = 2017; % Calendar end date (currently 2017 when GMF was sampled)
run_time = endtime;

pramvals.D = D;
pramvals.dz = dz;
pramvals.endtime = endtime;
pramvals.end_date = end_date;
%% Parameters
lasl = 0; % low asl = 0.002 
iasl = 0; % mid asl = 0.01
hasl = 0; % high asl = 0.03
lksl = 0; % low ksl = 0.01
iksl = 0; % mid ksl = 0.1
hksl = 0; % high ksl = 0.2

% if lasl == 1
%     a0_sl = 0.002;
% elseif iasl == 1;  
%     a0_sl = 0.01;  
% elseif hasl == 1;
%     a0_sl = 0.03;
% end 
% 
% if lksl == 1
%     ksl_surf = 0.01;
% elseif iksl == 1
%     ksl_surf = 0.1;
% elseif hksl == 1
%     ksl_surf = 0.2;
% end 

%0.2;%param(2);    % Slow pool consumption rate (1/yr)
%param(15);      % Zero-depth advection constant - slow (m/yr)

% load param % table storing parameter values fit from MCMC to GMF data
kr_surf = param(1);     % Rapid pool consumption rate (1/yr)
ksl_surf = param(2);    % Slow pool consumption rate (1/yr)
kst_surf = param(3);    % Stable pool consumption rate (1/yr)
er_surf = param(4);     % Rapid pool carbon use efficiency (CUE)
esl_surf = param(5);    % Slow pool CUE
est_surf = param(6);    % Stable pool CUE
tr_surf = param(7);     % Transformation constant of rapid C to slow C
tsl_surf = param(8);    % Transformation constant of slow C to stable C
Z = param(9);           % Root depth scaling constant (m)
tau = param(10);        % Rate of decrease of a0 with depth (1/m)
rootC_total = param(11); % Depth-integrated root C production (mol/m^2yr)
litterC_total = param(12); % Litter C input (mol/m^2yr)
bts = param(13);        % Depth scaling term for bioturbation (mol/m^2yr)
a0_r = param(14);       % Zero-depth advection constant - rapid (m/yr)
a0_sl = param(15);      % Zero-depth advection constant - slow (m/yr)
a0_st = param(16);      % Zero-depth advection constant - stable (m/yr)
Db = param(17);         % Bulk diffusion constant in m2/year 
% Db = 0; % No bioturbation (use when testing to increase run speed)
Dg = param(18).*ones(Nz,1); %Soil CO2 diffusion constant (m^2/yr) 
% Dg = zeroes(Nz,1); % Uncomment to remove CO2 diffusion
alpha_r = param(19);    % Fractionation during respiraiton (rapid) - Set to be from topsoil incubation experiments
alpha_sl = param(20);   % Fractionation during respiraiton (slow)
alpha_st=  param(21);   % Fractionation during respiraiton (stable)
%% Parameters run 
run_param_table(1) = kr_surf;
run_param_table(2) = ksl_surf;
run_param_table(3) = kst_surf;
run_param_table(4) = er_surf;
run_param_table(5) = esl_surf;
run_param_table(6) = est_surf;
run_param_table(7) = tr_surf;
run_param_table(8) = tsl_surf;
run_param_table(9) = Z;
run_param_table(10) = tau;
run_param_table(11) = rootC_total;
run_param_table(12) = litterC_total;
run_param_table(13) = bts;
run_param_table(14) = a0_r;
run_param_table(15) = a0_sl;
run_param_table(16) = a0_st;
run_param_table(17) = Db;
run_param_table(18) = Dg(1);
run_param_table(19) = alpha_r;
run_param_table(20) = alpha_sl;
run_param_table(21) = alpha_st;
% SUP = 0;
% [new_runtime,SpUp] = spinup(pramvals,run_param_table,SUP);

% endtime = 3000; % Duration of run in years
% end_date = 2017; % Calendar end date (currently 2017 when GMF was sampled)
% run_time = new_runtime;


%% Bulk C diffusion (bioturbation)
Df = Db.*exp(-z.*bts); % Scale bioturbation with depth
% Bioturbation is equal for each pool (below) 
Df_r = Df; % Rapid bioturbation
Df_sl = Df; % Slow bioturbation
Df_st = Df; % Stable bioturbation 
%% CO2 diffusion
% CO2 diffusion constants for 13C, 12C, and 14C 
D_13C =Dg/1.0044;
D_12C= Dg;
D_14C =Dg/1.00868; % Wang 1993 
poros = 0.49; % Porosity
%% Lag time for litter input 
t_max = 3;
lag_l = 0; % lag time currently set to 0 (litter C is immediately available for consumption)
lag_r = 0;
%% Alpha
alpha_r  = 0.9998400; % From topsoil incubation experiments
alpha_sl = 0.9999999; % Constant at ~1
alpha_st = 0.9999999; % Constant at ~1

% % Uncomment is alpha is the value from parameter table 
alpha_r = param(19);
alpha_sl = param(20);
alpha_st=  param(20);

kr = kr_surf.*ones(Nz,1);
ksl = ksl_surf.*ones(Nz,1);
%         ksl = ksl_surf.*exp(-(z./Z)); % Change slow decomp constant based on root C w/ depth     
kst = kst_surf.*ones(Nz,1);
%         kst = kst_surf.*exp(-(z./Z)); % Change stable decomp constant based on root C w/ depth  
er = er_surf.*ones(Nz,1);
esl = esl_surf.*ones(Nz,1);
est = est_surf.*ones(Nz,1);
tr = tr_surf.*ones(Nz,1);
tsl = tsl_surf.*ones(Nz,1);   
        
%% Advection terms 
az_r = a0_r*(1-(tau.*z)); % Vertical distribution of rapid C due to advection (m/year)
az_sl = a0_sl*(1-(tau.*z)); % Vertical distribution of slow C due to advection (m/year)
az_st = a0_st*(1-(tau.*z)); % Vertical distribution of stable C due to advection (m/year)

a0_max = max(abs([a0_r a0_sl a0_st])); % Maximum of all advection terms (for CFL/setting time step)
%% Set time step based on advection/diffusion terms
F = 0.5;
dt_a =F*dz/a0_max; % Time steps in units of years (CFL)
dt_d = F*((dz^2)./Df)/2;
dt_diff = ((0.01^2)./Dg)/2;
dt_diff = dt_diff(1);

% Change the time step depending on if diffusion is running
if Db > 0   
    dt_a = min([dt_a; dt_d]);
end 

t_total = dt_a:dt_a:endtime;

%% Vegetation 13C/12C ratio
% d13Cp = -26.6;  % d13C of the vegetation (from Torn et al 2002)
d13Cp = -28.6;
RC = (Rpdb*((d13Cp/1000)+1)).*ones(1,length(t_total)); % Ratio of the vegetation 
%% Radiocarbon 
lambda = 1/8267;%log(2)/5730; %half life of 14C in years 
decay_term = exp(-lambda*dt_a);
%% Roots
%If atmospheric CO2 concentrations and d13C values change through time, then these
%values need to be updated during each iteration in the time loop
rootC_r_total = .6*rootC_total;% sum of rapid + slow = 0.4%rootC_total * 0.67; % 67% of root C is rapid
rootC_sl_total = .4*rootC_total;%rootC_total * 0.33; % 33% of root C is slow
rootC_st_total = 0; % No root C is stable 

% % Uncomment to make C root input 0
% root_input_r = 0; 
% root_input_sl = 0; 
% root_input_st = 0; 

rootC_r =  dt_a*((rootC_r_total/Z)* exp(-(z./Z))); 
rootC_13r = rootC_r * (RC(1)/(1+RC(1)));
rootC_12r = rootC_r * (1/(1+RC(1)));
                   
rootC_sl = dt_a*((rootC_sl_total/Z)*exp(-(z./Z)));
rootC_13sl = rootC_sl * (RC(1)/(1+RC(1)));
rootC_12sl = rootC_sl * (1/(1+RC(1)));

rootC_st = dt_a*((rootC_st_total/Z)*exp(-(z./Z)));
rootC_13st = rootC_st * (RC(1)/(1+RC(1)));
rootC_12st = rootC_st* (1/(1+RC(1)));

%% Atmospheric conditions
d13Ca = d13Cfit (length(t_total),run_time,end_date)'; %atmospheric d13C going back to 1692
pCO2 = pCO2fit(length(t_total),run_time,end_date)'; %atmospheric pCO2 going back to 1694
C14_atm = C14fit(length(t_total),run_time,end_date)'; % 14C in fmc (data from Reimer 2013)
% C14_atm = zeros(length(t_total),1); % to make atm 14C 0 
pCO2= (flipud(pCO2))';
d13Ca = (flipud(d13Ca))';
C14_atm = (flipud(C14_atm))'; 

% Warnings if errors in interpolation of atmospheric data
if any(isnan(pCO2))
    disp('Error: bad interpolation of pCO2 data')
end 
if any(isnan(d13Ca))
    disp('Error: bad interpolation of d13C data')
end       
if any(isnan(C14_atm))
    disp('Error: bad interpolation of 14C data')
end 
%% Initialize the model
SUP = 0;
initialize_soilC_model_LRS_3
%% Schubert and Jahren d13C discrimination calcs
mtype = 1; % 0 - hold pCO2 constant; 1 - recorded atmospheric values

A = 28.26;
B = 0.22;
C = 23.9;
% d13Cai = d13Ca(1789);
d13Cai = d13Ca(length(t_total)-1); % Atmospheric d13CO2 values at end date
d13Cpi = -28.26; % From HR litter data in 2015 with same/similar forest as GMF
D13Ctp = (d13Cai - (d13Cpi))/(1+(d13Cpi)/1000); % Initial discrimination given modern pCO2 and d13Cveg (-27.5)
D13Ci = ((A*B)*(pCO2(length(t_total)-1)+C))./(A+B*(pCO2(length(t_total)-1)+C));
%% Time and depth loop
count = 0;
for time = 1:length(t_total)-1% time loop
% % Changing surface conditions of (atmosphere and vegetation)
    if mtype == 0
        d13Cp(time) = (d13Ca(time) - D13Ctp)./(1+D13Ctp/1000);
    else 
        D13C_Schubert =((A*B)*(pCO2(time)+C))./(A+B*(pCO2(time)+C)); % discrimination based on pCO2 (2015)
        DD13C = D13C_Schubert - D13Ci;
        D13C(time) = D13Ctp + DD13C;
        d13Cp(time) = (d13Ca(time) - D13C(time))./(1+D13C(time)/1000);
    end 

    Aatm(time) = (((C14_atm(time)/1000)+1)/((0.975^2)/((1+(d13Ca(time)/1000))^2)))*Aabs; % 14C Activity of the atmosphere
    D14Cveg(time) = 1000*(((Aatm(time)/Aabs)*((0.975^2)/((1+(d13Ca(time)/1000))^2)))-1); %D14C of vegetation
    Ratm(time) = (d13Ca(time)/1000+1)*Rpdb; % 12C/13C ratio of atmosphere
    RC(time) = d13Cp(time)/1000*Rpdb + Rpdb; % 12C/13C ratio of veg carbon
    
    molCO2_atm(time) = pCO2(time)/.022400/10^6; % mol CO2 atm
    mol12_CO2_atm(time) = molCO2_atm(time)/(Ratm(time)+1); % mol 12CO2 atm
    mol13_CO2_atm(time) = (Ratm(time)*molCO2_atm(time))/(Ratm(time)+1); % mol 13CO2 atm
    mol14_CO2_atm(time) = (((C14_atm(time)/1000)+1)/((0.975^2)/((1+(d13Ca(time)/1000))^2)))*Aabs*molCO2_atm(time); % mol 14CO2 atm
    
   for depth = 1:Nz-1 
        % No change in k, e, or t with depth
  
      
% % % Litter input
    if time*dt_a < lag_l % When time < lag, no litter is added yet. 
            litterC_r = 0;
            litterC_sl = 0;
            litterC_st = 0;
    elseif time*dt_a > lag_l 
            litterC_r = (litterC_total*0.6)*dt_a/0.01;
            litterC_sl = (litterC_total*0.4)*dt_a/0.01;
            litterC_st = 0;%(litterC_total*0.1)*dt_a/dz;
    end 
            litterC_13r = (RC(time-floor(lag_l/dt_a)).*litterC_r)./(1+RC(time-floor(lag_l/dt_a)));
            litterC_12r = litterC_r.*(1/(RC(time-floor(lag_l/dt_a))+1)) ;
            litterC_14r = (((D14Cveg(time-floor(lag_l/dt_a))/1000)+1)/((0.975^2)/((1+(d13Cp(time-floor(lag_l/dt_a))/1000))^2)))*Aabs*litterC_r;      

            litterC_13sl = (RC(time-floor(lag_l/dt_a)).*litterC_sl)./(1+RC(time-floor(lag_l/dt_a)));
            litterC_12sl =litterC_sl.*(1/(RC(time-floor(lag_l/dt_a))+1)) ;
            litterC_14sl = (((D14Cveg(time-floor(lag_l/dt_a))/1000)+1)/((0.975^2)/((1+(d13Cp(time-floor(lag_l/dt_a))/1000))^2)))*Aabs*litterC_sl;

            litterC_13st = (RC(time-floor(lag_l/dt_a)).*litterC_st)./(1+RC(time-floor(lag_l/dt_a)));
            litterC_12st =litterC_st.*(1/(RC(time-floor(lag_l/dt_a))+1)) ;
            litterC_14st = (((D14Cveg(time-floor(lag_l/dt_a))/1000)+1)/((0.975^2)/((1+(d13Cp(time-floor(lag_l/dt_a))/1000))^2)))*Aabs*litterC_st;
            
            rootC_13r(depth) = rootC_r(depth) * (RC(time)/(1+RC(time)));
            rootC_12r(depth) = rootC_r(depth) * (1/(1+RC(time)));
            rootC_14r(depth) = (((D14Cveg(time)/1000)+1)/((0.975^2)/((1+(d13Cp(time)/1000))^2)))*Aabs*rootC_r(depth);

            rootC_13sl(depth) = rootC_sl(depth) * (RC(time)/(1+RC(time)));
            rootC_12sl(depth) = rootC_sl(depth) * (1/(1+RC(time)));
            rootC_14sl(depth) = (((D14Cveg(time)/1000)+1)/((0.975^2)/((1+(d13Cp(time)/1000))^2)))*Aabs*rootC_sl(depth);

            rootC_13st(depth) = rootC_st(depth) * (RC(time)/(1+RC(time)));
            rootC_12st(depth) = rootC_st(depth) * (1/(1+RC(time)));
            rootC_14st(depth) = (((D14Cveg(time)/1000)+1)/((0.975^2)/((1+(d13Cp(time)/1000))^2)))*Aabs*rootC_st(depth);

    if depth ==1
%% At the surface %%
       %% Components for the time and depth loop    
        C_consumed_r(depth,time)    = dt_a*kr(depth).*Crz(depth,time);
        C_consumed_r_13(depth,time) = dt_a*kr(depth).*Crz_13(depth,time); 
        C_consumed_r_12(depth,time) = dt_a*kr(depth).*Crz_12(depth,time);
        C_consumed_r_14(depth,time) = dt_a*kr(depth).*Crz_14(depth,time); 
        
        C_consumed_sl(depth,time)    = dt_a*ksl(depth).*Cslz(depth,time);
        C_consumed_sl_13(depth,time) = dt_a*ksl(depth).*Cslz_13(depth,time);
        C_consumed_sl_12(depth,time) = dt_a*ksl(depth).*Cslz_12(depth,time);
        C_consumed_sl_14(depth,time) = dt_a*ksl(depth).*Cslz_14(depth,time);
        
        C_consumed_st(depth,time)    = dt_a*kst(depth).*Cstz(depth,time);
        C_consumed_st_13(depth,time) = dt_a*kst(depth).*Cstz_13(depth,time);
        C_consumed_st_12(depth,time) = dt_a*kst(depth).*Cstz_12(depth,time);
        C_consumed_st_14(depth,time) = dt_a*kst(depth).*Cstz_14(depth,time); 
        
        CO2_lost_r(depth,time)   = dt_a*kr(depth).*Crz(depth,time)*(1-er(depth));
        CO2_lost_13r(depth,time) = dt_a*kr(depth).*Crz(depth,time)*(1-er(depth))*CO2_13fac_r(depth);
        CO2_lost_12r(depth,time) = dt_a*kr(depth).*Crz(depth,time)*(1-er(depth))*CO2_12fac_r(depth);
        CO2_lost_14r  (depth,time) = dt_a*kr(depth).*Crz_14(depth,time)*(1-er(depth));
         
        CO2_lost_sl(depth,time)   = dt_a*ksl(depth).*Cslz(depth,time)*(1-esl(depth));
        CO2_lost_13sl(depth,time) = dt_a*ksl(depth).*Cslz(depth,time)*(1-esl(depth))*CO2_13fac_sl(depth);
        CO2_lost_12sl(depth,time) = dt_a*ksl(depth).*Cslz(depth,time)*(1-esl(depth))*CO2_12fac_sl(depth);
        CO2_lost_14sl (depth,time) = dt_a*ksl(depth).*Cslz_14(depth,time)*(1-esl(depth));
        
        CO2_lost_st(depth,time)   = dt_a*kst(depth).*Cstz(depth,time)*(1-est(depth));
        CO2_lost_13st(depth,time) = dt_a*kst(depth).*Cstz(depth,time)*(1-est(depth))*CO2_13fac_st(depth);
        CO2_lost_12st(depth,time) = dt_a*kst(depth).*Cstz(depth,time)*(1-est(depth))*CO2_12fac_st(depth);
        CO2_lost_14st (depth,time) = dt_a*kst(depth).*Cstz_14(depth,time)*(1-est(depth));
    
        CO2_lost_total(depth,time)     = (CO2_lost_r(depth,time) + CO2_lost_sl(depth,time) +CO2_lost_st(depth,time) + (rootC_r(depth)+rootC_sl(depth)))/dt_a;
        CO2_lost_total_13(depth,time)  = (CO2_lost_13r(depth,time) +CO2_lost_13sl(depth,time) +CO2_lost_13st(depth,time) + (rootC_13r(depth)+rootC_13sl(depth)))/dt_a;
        CO2_lost_total_12(depth,time)  = (CO2_lost_12r(depth,time) +CO2_lost_12sl(depth,time) +CO2_lost_12st(depth,time) + (rootC_12r(depth)+rootC_12sl(depth)))/dt_a;
        CO2_lost_total_14(depth,time)  = (CO2_lost_14r(depth,time) +CO2_lost_14sl(depth,time) +CO2_lost_14st(depth,time) + (rootC_14r(depth)+rootC_14sl(depth)))/dt_a;       
        
        growth_r(depth,time) = dt_a*kr(depth).*Crz(depth,time)*er(depth);
        transformed_C_r_to_sl(depth,time) = growth_r(depth,time) * tr(depth);
        
        growth_r_13(depth,time) = (C_consumed_r_13(depth,time) - CO2_lost_13r(depth,time)); %13C flux into microbial growth is equal to the total 13C flux into microbes minus the 13C flux to CO2.
        growth_r_12(depth,time) = (C_consumed_r_12(depth,time) - CO2_lost_12r(depth,time));
        growth_r_14(depth,time) =  (C_consumed_r_14(depth,time) - CO2_lost_14r(depth,time)); %growth_r(depth,time)*Aabs;%
    
        transformed_C_r_to_sl_13(depth,time) = growth_r_13(depth,time) * tr(depth);  %Multiply growth by tr to get 13C flux transferred to next pool
        transformed_C_r_to_sl_12(depth,time) = growth_r_12(depth,time) * tr(depth);
        transformed_C_r_to_sl_14(depth,time) = growth_r_14(depth,time) * tr(depth); %transformed_C_r_to_sl(depth,time)*Aabs;%
        
        growth_sl(depth,time)   =  dt_a*ksl(depth).*Cslz(depth,time)*esl(depth);
        transformed_C_sl_to_st(depth,time) = growth_sl(depth,time) * tsl(depth);
        
        growth_sl_13(depth,time) = (C_consumed_sl_13(depth,time) - CO2_lost_13sl(depth,time));  %13C flux into microbial growth is equal to the total 13C flux into microbes minus the 13C flux to CO2.
        growth_sl_12(depth,time) = (C_consumed_sl_12(depth,time) - CO2_lost_12sl(depth,time));
        growth_sl_14(depth,time) = (C_consumed_sl_14(depth,time) - CO2_lost_14sl(depth,time)); %growth_sl(depth,time)*Aabs;%
        
        transformed_C_sl_to_st_13(depth,time) = growth_sl_13(depth,time) * tsl(depth);  %Multiply growth by tr to get 13C flux transferred to next pool
        transformed_C_sl_to_st_12(depth,time) = growth_sl_12(depth,time) * tsl(depth);
        transformed_C_sl_to_st_14(depth,time) = growth_sl_14(depth,time) * tsl(depth);%transformed_C_sl_to_st(depth,time)*Aabs;%

% Calculating the C, 12C, and 13C concentration in each of the 3 pools
      Crz(depth,time+1)     = Crz(depth,time)  - (CO2_lost_r(depth,time) + transformed_C_r_to_sl(depth,time)) +  rootC_r(depth)+ litterC_r;
      Crz_13(depth,time+1)  = Crz_13(depth,time)  - (CO2_lost_13r(depth,time) + transformed_C_r_to_sl_13(depth,time)) + rootC_13r(depth) +litterC_13r ; 
      Crz_12(depth,time+1)  = Crz_12(depth,time)  - (CO2_lost_12r(depth,time)+ transformed_C_r_to_sl_12(depth,time)) + rootC_12r(depth) + litterC_12r ; 
      Crz_14(depth,time+1)  = Crz_14(depth,time)*decay_term - (CO2_lost_14r(depth,time) + transformed_C_r_to_sl_14(depth,time)) + rootC_14r(depth) + litterC_14r ;
    
      Cslz(depth,time+1)    =  Cslz(depth,time)  - (CO2_lost_sl(depth,time) + transformed_C_sl_to_st(depth,time)) + transformed_C_r_to_sl(depth,time) + rootC_sl(depth) +litterC_sl;
      Cslz_13(depth,time+1) = Cslz_13(depth,time) - (CO2_lost_13sl(depth,time)+ transformed_C_sl_to_st_13(depth,time)) + transformed_C_r_to_sl_13(depth,time) + rootC_13sl(depth)+litterC_13sl;
      Cslz_12(depth,time+1) = Cslz_12(depth,time) - (CO2_lost_12sl(depth,time)+ transformed_C_sl_to_st_12(depth,time)) + transformed_C_r_to_sl_12(depth,time) + rootC_12sl(depth)+litterC_12sl;
      Cslz_14(depth,time+1) = Cslz_14(depth,time)*decay_term - (CO2_lost_14sl(depth,time)+ transformed_C_sl_to_st_14(depth,time)) + transformed_C_r_to_sl_14(depth,time)*decay_term + rootC_14sl(depth) +litterC_14sl;
      
      Cstz(depth,time+1)    = Cstz(depth,time) - (CO2_lost_st(depth,time)) + transformed_C_sl_to_st(depth,time) + litterC_st;
      Cstz_13(depth,time+1) = Cstz_13(depth,time) - (CO2_lost_13st(depth,time)) + transformed_C_sl_to_st_13(depth,time) + litterC_13st;
      Cstz_12(depth,time+1) = Cstz_12(depth,time) - (CO2_lost_12st(depth,time)) + transformed_C_sl_to_st_12(depth,time)+ litterC_12st;
      Cstz_14(depth,time+1) = Cstz_14(depth,time)*decay_term - (CO2_lost_14st(depth,time)) + transformed_C_sl_to_st_14(depth,time)*decay_term + litterC_14st*decay_term;
      
    elseif depth >1  
%% Not at the surface!!!! %% 
 %% Components for the time and depth loop    
        C_consumed_r(depth,time)    = dt_a*kr(depth).*Crz(depth,time);
        C_consumed_r_13(depth,time) = dt_a*kr(depth).*Crz_13(depth,time); 
        C_consumed_r_12(depth,time) = dt_a*kr(depth).*Crz_12(depth,time);
        C_consumed_r_14(depth,time) = dt_a*kr(depth).*Crz_14(depth,time); %C_consumed_r(depth,time)*Aabs;% dt_a*kr(depth).*Crz(depth,time)*Aabs;%
        
        C_consumed_sl(depth,time)    = dt_a*ksl(depth).*Cslz(depth,time);
        C_consumed_sl_13(depth,time) = dt_a*ksl(depth).*Cslz_13(depth,time);
        C_consumed_sl_12(depth,time) = dt_a*ksl(depth).*Cslz_12(depth,time);
        C_consumed_sl_14(depth,time) = dt_a*ksl(depth).*Cslz_14(depth,time); %C_consumed_sl(depth,time)*Aabs;%dt_a*ksl(depth).*Cslz(depth,time)*Aabs;%
        
        C_consumed_st(depth,time)    = dt_a*kst(depth).*Cstz(depth,time);
        C_consumed_st_13(depth,time) = dt_a*kst(depth).*Cstz_13(depth,time);
        C_consumed_st_12(depth,time) = dt_a*kst(depth).*Cstz_12(depth,time);
        C_consumed_st_14(depth,time) = dt_a*kst(depth).*Cstz_14(depth,time); %C_consumed_st(depth,time)*Aabs;%dt_a*kst(depth).*Cstz(depth,time)*Aabs;%
        
        CO2_lost_r(depth,time)   = dt_a*kr(depth).*Crz(depth,time)*(1-er(depth));
        CO2_lost_13r(depth,time) = dt_a*kr(depth).*Crz(depth,time)*(1-er(depth))*CO2_13fac_r(depth);
        CO2_lost_12r(depth,time) = dt_a*kr(depth).*Crz(depth,time)*(1-er(depth))*CO2_12fac_r(depth);
        CO2_lost_14r  (depth,time) = dt_a*kr(depth).*Crz_14(depth,time)*(1-er(depth));%CO2_lost_r(depth,time) * Aabs;
         
        CO2_lost_sl(depth,time)   = dt_a*ksl(depth).*Cslz(depth,time)*(1-esl(depth));
        CO2_lost_13sl(depth,time) = dt_a*ksl(depth).*Cslz(depth,time)*(1-esl(depth))*CO2_13fac_sl(depth);
        CO2_lost_12sl(depth,time) = dt_a*ksl(depth).*Cslz(depth,time)*(1-esl(depth))*CO2_12fac_sl(depth);
        CO2_lost_14sl (depth,time) = dt_a*ksl(depth).*Cslz_14(depth,time)*(1-esl(depth)); %CO2_lost_sl(depth,time) *Aabs ;%
        
        CO2_lost_st(depth,time)   = dt_a*kst(depth).*Cstz(depth,time)*(1-est(depth));
        CO2_lost_13st(depth,time) = dt_a*kst(depth).*Cstz(depth,time)*(1-est(depth))*CO2_13fac_st(depth);
        CO2_lost_12st(depth,time) = dt_a*kst(depth).*Cstz(depth,time)*(1-est(depth))*CO2_12fac_st(depth);
        CO2_lost_14st (depth,time) = dt_a*kst(depth).*Cstz_14(depth,time)*(1-est(depth)); %CO2_lost_st(depth,time) *Aabs;%
         
        growth_r(depth,time) = dt_a*kr(depth).*Crz(depth,time)*er(depth);
        transformed_C_r_to_sl(depth,time) = growth_r(depth,time) * tr(depth);
        
        growth_r_13(depth,time) = (C_consumed_r_13(depth,time) - CO2_lost_13r(depth,time)); %13C flux into microbial growth is equal to the total 13C flux into microbes minus the 13C flux to CO2.
        growth_r_12(depth,time) = (C_consumed_r_12(depth,time) - CO2_lost_12r(depth,time));
        growth_r_14(depth,time) = (C_consumed_r_14(depth,time) - CO2_lost_14r(depth,time)); %growth_r(depth,time)*Aabs;%
    
        transformed_C_r_to_sl_13(depth,time) = growth_r_13(depth,time) * tr(depth);  %Multiply growth by tr to get 13C flux transferred to next pool
        transformed_C_r_to_sl_12(depth,time) = growth_r_12(depth,time) * tr(depth);
        transformed_C_r_to_sl_14(depth,time) = growth_r_14(depth,time) * tr(depth); %transformed_C_r_to_sl(depth,time)*Aabs;%
        
        growth_sl(depth,time)   =  dt_a*ksl(depth).*Cslz(depth,time)*esl(depth);
        transformed_C_sl_to_st(depth,time) = growth_sl(depth,time) * tsl(depth);
        
        growth_sl_13(depth,time) = (C_consumed_sl_13(depth,time) - CO2_lost_13sl(depth,time));  %13C flux into microbial growth is equal to the total 13C flux into microbes minus the 13C flux to CO2.
        growth_sl_12(depth,time) = (C_consumed_sl_12(depth,time) - CO2_lost_12sl(depth,time));
        growth_sl_14(depth,time) = (C_consumed_sl_14(depth,time) - CO2_lost_14sl(depth,time)); %growth_sl(depth,time)*Aabs;%
        
        transformed_C_sl_to_st_13(depth,time) = growth_sl_13(depth,time) * tsl(depth);  %Multiply growth by tr to get 13C flux transferred to next pool
        transformed_C_sl_to_st_12(depth,time) = growth_sl_12(depth,time) * tsl(depth);
        transformed_C_sl_to_st_14(depth,time) = growth_sl_14(depth,time) * tsl(depth); %transformed_C_sl_to_st(depth,time)*Aabs;%
        
% % %         Total CO2 lost  (CO2 production term) 
        CO2_lost_total(depth,time)     = (CO2_lost_r(depth,time) + CO2_lost_sl(depth,time) +CO2_lost_st(depth,time) + (rootC_r(depth)+rootC_sl(depth)))/dt_a;
        CO2_lost_total_13(depth,time)  = (CO2_lost_13r(depth,time) +CO2_lost_13sl(depth,time) +CO2_lost_13st(depth,time) + (rootC_13r(depth)+rootC_13sl(depth)))/dt_a;
        CO2_lost_total_12(depth,time)  = (CO2_lost_12r(depth,time) +CO2_lost_12sl(depth,time) +CO2_lost_12st(depth,time) + (rootC_12r(depth)+rootC_12sl(depth)))/dt_a;
        CO2_lost_total_14(depth,time)  = (CO2_lost_14r(depth,time) +CO2_lost_14sl(depth,time) +CO2_lost_14st(depth,time) + (rootC_14r(depth)+rootC_14sl(depth)))/dt_a;
        
% % %         Advective gradient coming down from the above depth    
        Advgradient_13r(depth,time)  =  az_r(depth)*dt_a.*(Crz_13(depth,time)-Crz_13(depth-1,time))/dz;
        Advgradient_r(depth,time)  =  az_r(depth)*dt_a.*(Crz(depth,time)-Crz(depth-1,time))/dz;
        Advgradient_12r(depth,time)  =  az_r(depth)*dt_a.*(Crz_12(depth,time)-Crz_12(depth-1,time))/dz;
        Advgradient_14r(depth,time)  =  az_r(depth)*dt_a.*(Crz_14(depth,time)-Crz_14(depth-1,time))/dz; %Advgradient_r(depth,time)*Aabs;%

        Advgradient_13sl(depth,time)  = az_sl(depth)*dt_a.*(Cslz_13(depth,time)-Cslz_13(depth-1,time))/dz;
        Advgradient_sl(depth,time)  =  az_sl(depth)*dt_a.*(Cslz(depth,time)-Cslz(depth-1,time))/dz;
        Advgradient_12sl(depth,time) =  az_sl(depth)*dt_a.*(Cslz_12(depth,time)-Cslz_12(depth-1,time))/dz;
        Advgradient_14sl(depth,time) =  az_sl(depth)*dt_a.*(Cslz_14(depth,time)-Cslz_14(depth-1,time))/dz; %Advgradient_sl(depth,time)*Aabs;%

        Advgradient_13st(depth,time)  =  az_st(depth)*dt_a.*(Cstz_13(depth,time)-Cstz_13(depth-1,time))/dz;
        Advgradient_st(depth,time)  =  az_st(depth)*dt_a.*(Cstz(depth,time)-Cstz(depth-1,time))/dz;
        Advgradient_12st(depth,time) =  az_st(depth)*dt_a.*(Cstz_12(depth,time)-Cstz_12(depth-1,time))/dz;
        Advgradient_14st(depth,time) =  az_st(depth)*dt_a.*(Cstz_14(depth,time)-Cstz_14(depth-1,time))/dz;% Advgradient_st(depth,time)*Aabs;%

% % %               Bioturbation/diffusion
        Cdiff_r(depth,time)   = dt_a*((Crz(depth-1,time)    -   2*Crz(depth,time) + Crz(depth+1,time))/(dz^2))*Df_r(depth);
        Cdiff_13r(depth,time) = dt_a*((Crz_13(depth-1,time) - 2*Crz_13(depth,time) + Crz_13(depth+1,time))/(dz^2))*Df_r(depth);
        Cdiff_12r(depth,time) = dt_a*((Crz_12(depth-1,time) - 2*Crz_12(depth,time) + Crz_12(depth+1,time))/(dz^2))*Df_r(depth);
        Cdiff_14r(depth,time) = dt_a*((Crz_14(depth-1,time) - 2*Crz_14(depth,time) + Crz_14(depth+1,time))/(dz^2))*Df_r(depth);
        
        Cdiff_sl(depth,time)   = dt_a*((Cslz(depth-1,time)    - 2*Cslz(depth,time) + Cslz(depth+1,time))/(dz^2))*Df_sl(depth);
        Cdiff_13sl(depth,time) = dt_a*((Cslz_13(depth-1,time) - 2*Cslz_13(depth,time) + Cslz_13(depth+1,time))/(dz^2))*Df_sl(depth);
        Cdiff_12sl(depth,time) = dt_a*((Cslz_12(depth-1,time) - 2*Cslz_12(depth,time) + Cslz_12(depth+1,time))/(dz^2))*Df_sl(depth);
        Cdiff_14sl(depth,time) = dt_a*((Cslz_14(depth-1,time) - 2*Cslz_14(depth,time) + Cslz_14(depth+1,time))/(1*dz^2))*Df_sl(depth);
        
        Cdiff_st(depth,time)   = dt_a*((Cstz(depth-1,time)    - 2*Cstz(depth,time) + Cstz(depth+1,time))/(dz^2))*Df_st(depth);
        Cdiff_13st(depth,time) = dt_a*((Cstz_13(depth-1,time) - 2*Cstz_13(depth,time) + Cstz_13(depth+1,time))/(dz^2))*Df_st(depth);
        Cdiff_12st(depth,time) = dt_a*((Cstz_12(depth-1,time) - 2*Cstz_12(depth,time) + Cstz_12(depth+1,time))/(dz^2))*Df_st(depth);
        Cdiff_14st(depth,time) = dt_a*((Cstz_14(depth-1,time) - 2*Cstz_14(depth,time) + Cstz_14(depth+1,time))/(dz^2))*Df_st(depth);
       

      Crz(depth,time+1)     = Crz(depth,time)  + Cdiff_r(depth,time) - Advgradient_r(depth,time) - (CO2_lost_r(depth,time) + transformed_C_r_to_sl(depth,time)) +  rootC_r(depth); 
      Crz_13(depth,time+1)  = Crz_13(depth,time) + Cdiff_13r(depth,time) - Advgradient_13r(depth,time) - (CO2_lost_13r(depth,time) + transformed_C_r_to_sl_13(depth,time)) +rootC_13r(depth) ; 
      Crz_12(depth,time+1)  = Crz_12(depth,time) + Cdiff_12r(depth,time) - Advgradient_12r(depth,time) - (CO2_lost_12r(depth,time)+ transformed_C_r_to_sl_12(depth,time)) +rootC_12r(depth) ; 
      Crz_14(depth,time+1)  = Crz_14(depth,time)*decay_term + Cdiff_14r(depth,time)*decay_term - Advgradient_14r(depth,time)*decay_term - (CO2_lost_14r(depth,time)+ transformed_C_r_to_sl_14(depth,time)) +rootC_14r(depth); 
      
      Cslz(depth,time+1)    = Cslz(depth,time) + Cdiff_sl(depth,time) - Advgradient_sl(depth,time) -(CO2_lost_sl(depth,time) + transformed_C_sl_to_st(depth,time)) + transformed_C_r_to_sl(depth,time) + rootC_sl(depth);
      Cslz_13(depth,time+1) = Cslz_13(depth,time) + Cdiff_13sl(depth,time) - Advgradient_13sl(depth,time) -(CO2_lost_13sl(depth,time)+ transformed_C_sl_to_st_13(depth,time)) + transformed_C_r_to_sl_13(depth,time) + rootC_13sl(depth);
      Cslz_12(depth,time+1) = Cslz_12(depth,time) + Cdiff_12sl(depth,time) - Advgradient_12sl(depth,time) -(CO2_lost_12sl(depth,time)+ transformed_C_sl_to_st_12(depth,time)) + transformed_C_r_to_sl_12(depth,time) + rootC_12sl(depth);
      Cslz_14(depth,time+1) = Cslz_14(depth,time)*decay_term + Cdiff_14sl(depth,time)*decay_term - Advgradient_14sl(depth,time)*decay_term -(CO2_lost_14sl(depth,time)+ transformed_C_sl_to_st_14(depth,time)) + transformed_C_r_to_sl_14(depth,time)*decay_term + rootC_14sl(depth);
      
      Cstz(depth,time+1)    = Cstz(depth,time) +    Cdiff_st(depth,time) - Advgradient_st(depth,time) - (CO2_lost_st(depth,time)) + transformed_C_sl_to_st(depth,time);
      Cstz_13(depth,time+1) = Cstz_13(depth,time) + Cdiff_13st(depth,time) - Advgradient_13st(depth,time) -(CO2_lost_13st(depth,time)) + transformed_C_sl_to_st_13(depth,time);
      Cstz_12(depth,time+1) = Cstz_12(depth,time) + Cdiff_12st(depth,time) - Advgradient_12st(depth,time) -(CO2_lost_12st(depth,time)) + transformed_C_sl_to_st_12(depth,time);
      Cstz_14(depth,time+1) = Cstz_14(depth,time)*decay_term + Cdiff_14st(depth,time)*decay_term - Advgradient_14st(depth,time)*decay_term -(CO2_lost_14st(depth,time)) + transformed_C_sl_to_st_14(depth,time)*decay_term;
    end   
    
% % %         Calculating the ratio of each pool
        RCrz(depth,time+1)   = Crz_13(depth,time+1)./Crz_12(depth,time+1);
        RCslz(depth,time+1)  = Cslz_13(depth,time+1)./Cslz_12(depth,time+1);
        RCstz(depth,time+1)  = Cstz_13(depth,time+1)./Cstz_12(depth,time+1);        
% % %         Calculating d13C of soil
        d13Crz(depth,time+1)  = ((RCrz(depth,time+1)./Rpdb) -1)*1000;
        d13Cslz(depth,time+1)  = ((RCslz(depth,time+1)./Rpdb) -1)*1000;
        d13Cstz(depth,time+1)  = ((RCstz(depth,time+1)./Rpdb) -1)*1000;  
% % %         Calculating Delta 14C of soil        
        DC14_r(depth,time+1) =  ((((Crz_14(depth,time+1) ./(Crz(depth,time+1)))/(Aabs)).*((0.975^2)./((1+(d13Crz(depth,time+1)./1000)).^2)))-1).*1000;
        DC14_sl(depth,time+1) = ((((Cslz_14(depth,time+1) ./(Cslz(depth,time+1)))/(Aabs)).*((0.975^2)./((1+(d13Cslz(depth,time+1)./1000)).^2)))-1).*1000;
        DC14_st(depth,time+1) = ((((Cstz_14(depth,time+1) ./(Cstz(depth,time+1)))/(Aabs)).*((0.975^2)./((1+(d13Cstz(depth,time+1)./1000)).^2)))-1).*1000;
        
         if depth == Nz-1  % This is a lower boundary condition/no flow condition 
            Crz(depth+1,time+1) = Crz(depth,time+1);
            Cslz(depth+1,time+1) = Cslz(depth,time+1);
            Cstz(depth+1,time+1) = Cstz(depth,time+1);
            Crz_13(depth+1,time+1) = Crz_13(depth,time+1);
            Cslz_13(depth+1,time+1) = Cslz_13(depth,time+1);
            Cstz_13(depth+1,time+1) = Cstz_13(depth,time+1);
            Crz_12(depth+1,time+1) = Crz_12(depth,time+1);
            Cslz_12(depth+1,time+1) = Cslz_12(depth,time+1);
            Cstz_12(depth+1,time+1) = Cstz_12(depth,time+1);
            Crz_14(depth+1,time+1) = Crz_14(depth,time+1);
            Cslz_14(depth+1,time+1) = Cslz_14(depth,time+1);
            Cstz_14(depth+1,time+1) = Cstz_14(depth,time+1);
            RCrz(depth+1,time+1) = RCrz(depth,time+1);
            RCslz(depth+1,time+1) = RCslz(depth,time+1);
            RCstz(depth+1,time+1) = RCstz(depth,time+1);
            d13Crz(depth+1,time+1)  =  d13Crz(depth,time+1);
            d13Cslz(depth+1,time+1) =  d13Cslz(depth,time+1);
            d13Cstz(depth+1,time+1) =  d13Cstz(depth,time+1); 
            DC14_r(depth+1,time+1) = DC14_r(depth,time+1);
            DC14_sl(depth+1,time+1) = DC14_sl(depth,time+1);
            DC14_st(depth+1,time+1) = DC14_st(depth,time+1);
         end
   end
% % %     Update mass balance terms
    
    CO2_12fac_r = 1./(1 + (alpha_r*RCrz(:,time+1)));
    CO2_13fac_r = alpha_r*RCrz(:,time+1)./(1+alpha_r*RCrz(:,time+1));
    CO2_12fac_sl = 1./(1 + (alpha_sl*RCslz(:,time+1)));
    CO2_13fac_sl = alpha_sl*RCslz(:,time+1)./(1+alpha_sl*RCslz(:,time+1));
    CO2_12fac_st = 1./(1 + (alpha_st*RCstz(:,time+1)));
    CO2_13fac_st = alpha_st*RCstz(:,time+1)./(1+alpha_st*RCstz(:,time+1));

    % Find steady state CO2 profile after each time step
% [CO2_old,CO2_old_12,CO2_old_13,CO2_old_14] = SS_CO2(time,dz,Nz,Dg,poros,dt_diff, Aabs, d13Ca, pCO2, C14_atm, ...
%                                                     Ratm,CO2_lost_total,CO2_lost_total_12,CO2_lost_total_13,CO2_lost_total_14,...
%                                                     CO2_conc,CO2_conc_12,CO2_conc_13,CO2_conc_14);
%     % Update CO2 concentration profiles
%     CO2_conc(:,time) = CO2_old;
%     CO2_conc_12(:,time) = CO2_old_12;
%     CO2_conc_13(:,time) = CO2_old_13;
%     CO2_conc_14(:,time) = CO2_old_14;
    C_total(:,time+1) = Crz(:,time+1) + Cslz(:,time+1) + Cstz(:,time+1);
    FC_rapid(:,time+1)  = Crz(:,time+1)./C_total(:,time+1);
    FC_slow(:,time+1)  = Cslz(:,time+1)./C_total(:,time+1);
    FC_stable(:,time+1)= Cstz(:,time+1)./C_total(:,time+1);
    
    C14_bulk(:,time+1) = (DC14_r(:,time+1) .* FC_rapid(:,time+1) )+ ... 
                    (DC14_sl(:,time+1) .* FC_slow(:,time+1) )+...
                    (DC14_st(:,time+1).* FC_stable(:,time+1));
     C14_diff_max(time+1) = norm(C14_bulk(:,time+1) - C14_bulk(:,time))/norm(C14_bulk(:,time+1));           
     count = count+1;
end 
%% Post loop calculations
Crz_14_end =  DC14_r(:,end);
Cslz_14_end = DC14_sl(:,end);
Cstz_14_end = DC14_st(:,end);

Crz_end  =  Crz_13(:,end)  + Crz_12(:,end);
Cslz_end =  Cslz_13(:,end) + Cslz_12(:,end);
Cstz_end =  Cstz_13(:,end) + Cstz_12(:,end);

C_total_end_i = Crz_end + Cslz_end  + Cstz_end;

FC_rapid_end  = Crz_end./C_total_end_i;
FC_slow_end   = Cslz_end./C_total_end_i;
FC_stable_end = Cstz_end./C_total_end_i;

d13C_bulk_end_i   = (d13Crz(:,end).*FC_rapid_end) + ...
                    (d13Cslz(:,end).*FC_slow_end) + ...
                    (d13Cstz(:,end).*FC_stable_end);

C14_bulk_end_i =    (Crz_14_end .* FC_rapid_end )+ ... 
                    (Cslz_14_end .* FC_slow_end )+...
                    (Cstz_14_end.* FC_stable_end);
                
CO2_total_end  = CO2_conc_12(:,time) + CO2_conc_13(:,time); % CO2 in the profile at the end (in mol/m3)
CO2_total_end_i = CO2_total_end.*0.022400*10^6; % CO2 at the end in ppmv 
d13C_CO2_end_i = ((CO2_conc_13(:,time)./CO2_conc_12(:,time))./Rpdb-1)*1000;
%% Respired CO2 profiles
RRpro = sum((CO2_lost_total(:,time).*dz)) ; %CO2 flux/Resp rate in mol/m2year
CO2_lost_end_r = CO2_lost_r(:,time);
CO2_lost_end_sl = CO2_lost_sl(:,time);
CO2_lost_end_st = CO2_lost_st(:,time);
CO2_lost_end_tot = CO2_lost_end_r + CO2_lost_end_sl  + CO2_lost_end_st;

CO2_lost_total_13_end = CO2_lost_13r(:,time) +CO2_lost_13sl(:,time) +CO2_lost_13st(:,time);
CO2_lost_total_12_end = CO2_lost_12r(:,time) +CO2_lost_12sl(:,time) +CO2_lost_12st(:,time);
CO2_lost_total_14_end = CO2_lost_14r(:,time) +CO2_lost_14sl(:,time) +CO2_lost_14st(:,time);

% % save resp profiles
resp_profiles = zeros(Nz,7);
resp_profiles(:,1) = CO2_lost_end_tot;
resp_profiles(:,2) = CO2_lost_end_r;
resp_profiles(:,3) = CO2_lost_end_sl;
resp_profiles(:,4) = CO2_lost_end_st;
resp_profiles(:,5) = CO2_lost_total_13_end;
resp_profiles(:,6) = CO2_lost_total_12_end;
resp_profiles(:,7) = CO2_lost_total_14_end;

%%
lstyle = '-';
% lcolor = 'k';
% figure(2)
%   subplot(1,2,1)
%     hold on
%     plot(C14_bulk_end_i,z,'-','Color',lcolor,'linewidth',1.5) ;
%     plot(Cslz_14_end,z,'--','Color',lcolor,'linewidth',1.5) ;
%     plot(Cstz_14_end,z,'-.','Color',lcolor,'linewidth',1.5) ;
%     set(gca, 'YDir', 'reverse','linewidth',1.5) ;
%     ylabel('Depth (m)','fontsize',12)
%     xlabel(['{\Delta}^{14}C (', char(8240), ')'],'fontsize',12,'FontWeight','bold')
%   subplot(1,2,2)
% % figure(2)
%     hold on
%     plot(C_total_end_i,z,lstyle,'Color',lcolor,'linewidth',1.25)
%     set(gca, 'YDir', 'reverse','linewidth',1.5) ;
%     ylabel('Depth (m)','fontsize',12)
%     xlabel('SOC (mol/m^{3})','fontsize',12,'FontWeight', 'bold')
%%

% 
% if lasl == 1
%     lstyle = '--';
%     atit = 'Low a_{sl}';
%     atit2 = 'low_asl';
% elseif iasl == 1;  
%     lstyle = '-.';
%     atit = 'Mid a_{sl}';
%     atit2 = 'mid_asl';    
% elseif hasl == 1;
%     lstyle = '-';
%     atit = 'High a_{sl}';
%     atit2 = 'high_asl';
% end 
% 
% if lksl == 1
%     lcolor = [0.2 0.6 1.0];
%     btit = 'Low k_{sl}';
%     btit2 = 'low_ksl';
% elseif iksl == 1
%     lcolor = [0.0 0.6 0.6];
%     btit = 'Mid k_{sl}';
%     btit2 = 'mid_ksl';
% elseif hksl == 1
%     lcolor = [0.2 0.6 0.0];
%     btit = 'High k_{sl}';
%     btit2 = 'high_ksl';
% end 
% tit = sprintf('%s and %s',atit,btit);
% fname1 = sprintf('pools_%s_%s.fig',atit,btit);
% fname2 = sprintf('%s_%s.fig',atit,btit);
% % fname3 = sprintf('%s_%s.fig',atit,btit);
% figure
% suptitle(tit)
% %   subplot(1,2,2)
%     hold on
%     plot(Crz_end,z,lstyle,'Color','r','linewidth',1.5) ;
%     plot(Cslz_end,z,lstyle,'Color','b','linewidth',1.5) ;
%     plot(Cstz_end,z,lstyle,'Color','g','linewidth',1.5) ;
%     set(gca, 'YDir', 'reverse','linewidth',1.5) ;
%     legend('Rapid','Slow','Stable')
%     ylabel('Depth (m)','fontsize',12)
%     xlabel('SOC (mol/m^{3})','fontsize',12,'FontWeight', 'bold')
% 
% saveas(gca,fname1);
% 
% figure
% suptitle(tit)
%   subplot(1,4,1)
%     hold on
%     plot(C14_bulk_end_i,z,'-','Color',lcolor,'linewidth',1.5) ;
%     set(gca, 'YDir', 'reverse','linewidth',1.5) ;
%     xlabel(['Bulk {\Delta}^{14}C (', char(8240), ')'],...
%         'fontsize',12,'FontWeight','bold')
%     ylabel('Depth (m)','fontsize',12,'FontWeight','bold')
%   subplot(1,4,2)
%     hold on
%     plot(Crz_14_end,z,'-','Color',lcolor,'linewidth',1.5) ;
%     set(gca, 'YDir', 'reverse','linewidth',1.5) ;    
%     xlabel(['Rapid {\Delta}^{14}C (', char(8240), ')'],...
%         'fontsize',12,'FontWeight','bold')
%   subplot(1,4,3)
%     hold on    
%     plot(Cslz_14_end,z,'-','Color',lcolor,'linewidth',1.5) ;
%     set(gca, 'YDir', 'reverse','linewidth',1.5) ;
%     xlabel(['Slow {\Delta}^{14}C (', char(8240), ')'],...
%         'fontsize',12,'FontWeight','bold')
%    subplot(1,4,4)
%     hold on   
%     plot(Cstz_14_end,z,'-','Color',lcolor,'linewidth',1.5) ;
%     set(gca, 'YDir', 'reverse','linewidth',1.5) ;
% 
%     xlabel(['Stable {\Delta}^{14}C (', char(8240), ')'],...
%         'fontsize',12,'FontWeight','bold')
% % if hksl == 1;
% %     legend('Low k_{sl}', 'Mid k_{sl}', 'High k_{sl}')
% % end 
% saveas(gca,fname2);
% 

% if hasl == 1
%     figure
%     suptitle('High a_{sl} values')
%       subplot(1,4,1)
%         hold on
%         plot(C14_bulk_end_i,z,'-','Color',lcolor,'linewidth',1.5) ;
%         set(gca, 'YDir', 'reverse','linewidth',1.5) ;
%         xlabel(['Bulk {\Delta}^{14}C (', char(8240), ')'],...
%             'fontsize',12,'FontWeight','bold')
%         ylabel('Depth (m)','fontsize',12,'FontWeight','bold')
%       subplot(1,4,2)
%         hold on
%         plot(Crz_14_end,z,'-','Color',lcolor,'linewidth',1.5) ;
%         set(gca, 'YDir', 'reverse','linewidth',1.5) ;    
%         xlabel(['Rapid {\Delta}^{14}C (', char(8240), ')'],...
%             'fontsize',12,'FontWeight','bold')
%       subplot(1,4,3)
%         hold on    
%         plot(Cslz_14_end,z,'-','Color',lcolor,'linewidth',1.5) ;
%         set(gca, 'YDir', 'reverse','linewidth',1.5) ;
%         xlabel(['Slow {\Delta}^{14}C (', char(8240), ')'],...
%             'fontsize',12,'FontWeight','bold')
%        subplot(1,4,4)
%         hold on   
%         plot(Cstz_14_end,z,'-','Color',lcolor,'linewidth',1.5) ;
%         set(gca, 'YDir', 'reverse','linewidth',1.5) ;
% 
%         xlabel(['Stable {\Delta}^{14}C (', char(8240), ')'],...
%             'fontsize',12,'FontWeight','bold')
%     if hksl == 1;
%         legend('Low k_{sl}', 'Mid k_{sl}', 'High k_{sl}')
%     end 
%     saveas(gca,fname2);
% end 

%%
% figure(3)
% plot(C14_bulk_end_i,z)
% set(gca,'YDir','reverse')
% hold on
% figure(4)
% plot(C_total_end_i,z)
% set(gca,'YDir','reverse')
% hold on
%%
lcolor = 'r';
lstyle = '-.';
    figure(1)
    % subplot(1,3,3)
        hold on
        plot(d13C_bulk_end_i,z,lstyle,'Color',lcolor,'linewidth',1.5) ;
        set(gca, 'YDir', 'reverse','linewidth',1.5) ;
        ylabel('Depth (m)','fontsize',12)
        xlabel(['\delta^{13}C (', char(8240), ')'],'fontsize',12,'FontWeight', 'bold')
        