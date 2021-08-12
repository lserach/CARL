function [SS_runtime,SpUp] = spinup(pramvals,run_param_table,SUP)
% Find steady state model duration 
%% Model profile depth parameters
D = pramvals.D  ; % use for when only looking at top 1 meter 
dz = pramvals.dz; %The size of each cell in m (0.01 m = 1 cm)
z = (0:dz:D)';  % Depth vector 
Nz = D/dz+1; %Number of nodes in the z direction (depth)
Aabs = 10e-12; 
An = 6.0221409e23; % Avogadros number
Rpdb = 0.0112372; % In per mil, PDB standard ratio of 13C/12C

endtime = pramvals.endtime;
end_date = pramvals.end_date;
run_time = endtime;

param = run_param_table;
%% Parameters
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
Dg = param(18).*ones(Nz,1); %Soil CO2 diffusion constant (m^2/yr) 
alpha_r = param(19);    % Fractionation during respiraiton (rapid) - Set to be from topsoil incubation experiments
alpha_sl = param(20);   % Fractionation during respiraiton (slow)
alpha_st=  param(21);   % Fractionation during respiraiton (stable)
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
% alpha_r  = 0.9998400; % From topsoil incubation experiments
% alpha_sl = 0.9999999; % Constant at ~1
% alpha_st = 0.9999999; % Constant at ~1

% % Uncomment is alpha is the value from parameter table 
% alpha_r = param(19);
% alpha_sl = param(20);
% alpha_st=  param(21);

%% Advection terms 
az_r = a0_r*(1-(tau.*z)); % Vertical distribution of rapid C due to advection (m/year)
az_sl = a0_sl*(1-(tau.*z)); % Vertical distribution of slow C due to advection (m/year)
az_st = a0_st*(1-(tau.*z)); % Vertical distribution of stable C due to advection (m/year)

a0_max = max(abs([a0_r a0_sl a0_st])); % Maximum of all advection terms (for CFL/setting time step)
%% Set time step based on advection/diffusion terms
F = 0.1;
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
d13Cp = -26.6;  % d13C of the vegetation (from Torn et al 2002)
% d13Cp = -27.5;
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
C14_atm = C14fit(length(t_total),run_time,end_date)'; % 14C in D14C (data from Reimer 2013)
C14_atm = zeros(length(t_total),1);
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
initialize_soilC_model_LRS_3
%% Schubert and Jahren d13C discrimination calcs
mtype = 0; % 0 - hold pCO2 constant; 1 - recorded atmospheric values

A = 28.26;
B = 0.22;
C = 23.9;
% d13Cai = d13Ca(1789);
d13Cai = d13Ca(length(t_total)-1); % Atmospheric d13CO2 values at end date
d13Cpi = -28.26; % From HR litter data in 2015 with same/similar forest as GMF
D13Ctp = (d13Cai - (d13Cpi))/(1+(d13Cpi)/1000); % Initial discrimination given modern pCO2 and d13Cveg (-27.5)
D13Ci = ((A*B)*(pCO2(length(t_total)-1)+C))./(A+B*(pCO2(length(t_total)-1)+C));
        % No change in k, e, or t with depth
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
        
%% Time and depth loop
count = 0;
SS_check = 1;
while SS_check > 0.5 
for time = 1:length(t_total)-1% time loop
   if SS_check ==0
       break
   end 
% % Changing surface conditions of (atmosphere and vegetation)
    if mtype == 0
        d13Cp(time) = (d13Ca(time) - D13Ctp)./(1+D13Ctp/1000);
        d13Cp(time) = -27.5;
    else 
        D13C_Schubert =((A*B)*(pCO2(time)+C))./(A+B*(pCO2(time)+C)); % discrimination based on pCO2 (2015)
        DD13C = D13C_Schubert - D13Ci;
        D13C(time) = D13Ctp + DD13C;
        d13Cp(time) = (d13Ca(time) - D13C(time))./(1+D13C(time)/1000);
    end 

    Aatm(time) = (((C14_atm(time)/1000)+1)/((0.975^2)/((1+(d13Ca(time)/1000))^2)))*Aabs; % 14C Activity of the atmosphere
    D14Cveg(time)   = 1000*(((Aatm(time)/Aabs)*((0.975^2)/((1+(d13Ca(time)/1000))^2)))-1); %D14C of vegetation
%     D14Cveg(time) = 0;
    Ratm(time) = (d13Ca(time)/1000+1)*Rpdb; % 12C/13C ratio of atmosphere
    RC(time) = d13Cp(time)/1000*Rpdb + Rpdb; % 12C/13C ratio of veg carbon
    
    molCO2_atm(time) = pCO2(time)/.022400/10^6; % mol CO2 atm
    mol12_CO2_atm(time) = molCO2_atm(time)/(Ratm(time)+1); % mol 12CO2 atm
    mol13_CO2_atm(time) = (Ratm(time)*molCO2_atm(time))/(Ratm(time)+1); % mol 13CO2 atm
    mol14_CO2_atm(time) = (((C14_atm(time)/1000)+1)/((0.975^2)/((1+(d13Ca(time)/1000))^2)))*Aabs*molCO2_atm(time); % mol 14CO2 atm
    
   for depth = 1:Nz-1 
  
      
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
    
    d13C_bulk(:,time+1) = (d13Crz(:,time+1).*FC_rapid(:,time+1)) + ...
                    (d13Cslz(:,time+1).*FC_slow(:,time+1)) + ...
                    (d13Cstz(:,time+1).*FC_stable(:,time+1));
    C14_bulk(:,time+1) = (DC14_r(:,time+1) .* FC_rapid(:,time+1) )+ ... 
                    (DC14_sl(:,time+1) .* FC_slow(:,time+1) )+...
                    (DC14_st(:,time+1).* FC_stable(:,time+1));
    % Determine if at SS
    C14_diff_max(time+1) = norm(C14_bulk(:,time+1) - C14_bulk(:,time))/norm(C14_bulk(:,time+1));
    d13C_diff_max(time+1) = norm(d13C_bulk(:,time+1) - d13C_bulk(:,time))/norm(d13C_bulk(:,time+1));
    C_diff_max(time+1) = norm(C_total(:,time+1) - C_total(:,time))/norm(C_total(:,time+1));
    if C14_diff_max(time+1) <  1.0e-05
            if C_diff_max(time+1) < 1.0e-05
                SS_check = 0;
            end 

    end 
    end               
     count = count+1;
end 

SpUp.Crz_i  = Crz(:,end);
SpUp.Cslz_i = Cslz(:,end);
SpUp.Cstz_i = Cstz(:,end);

SpUp.Crz_12_i  = Crz_12(:,end);
SpUp.Cslz_12_i = Cslz_12(:,end);
SpUp.Cstz_12_i = Cstz_12(:,end);

SpUp.Crz_13_i  = Crz_13(:,end);
SpUp.Cslz_13_i = Cslz_13(:,end);
SpUp.Cstz_13_i = Cstz_13(:,end);

SpUp.Crz_14_i  = Crz_14(:,end);
SpUp.Cslz_14_i = Cslz_14(:,end);
SpUp.Cstz_14_i = Cstz_14(:,end);

SpUp.DC14_r_i = DC14_r(:,end);
SpUp.DC14_sl_i = DC14_sl(:,end);
SpUp.DC14_st_i = DC14_st(:,end);

SS_runtime = floor((time+1)*dt_a)+1;
