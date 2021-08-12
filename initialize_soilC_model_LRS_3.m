Rpdb = 0.0112372;
RC = (Rpdb*((d13Cp/1000)+1)).*ones(1,length(t_total)); %Ratio of the vegetation *ones(Nz,1)
%RC  = RCp*ones(length(t_total),1); % Initial ratio at the surface is the ratio of the vegetation
% d13Cp  = d13Cp.*ones(1,length(t_total));
% D13C_Schubert = zeros(1,length(t_total));
% d13Cp = zeros(1,length(t_total)); 
Aatm = zeros(1,length(t_total));
Ratm = zeros(1,length(t_total));
D14Cveg = zeros(1,length(t_total));
D13C = zeros(1,length(t_total));
% 
molCO2_atm = zeros(1,length(t_total));
mol12_CO2_atm = zeros(1,length(t_total));
mol13_CO2_atm = zeros(1,length(t_total));
mol14_CO2_atm = zeros(1,length(t_total));
%     
rootC_14r = zeros(Nz,length(t_total));
rootC_14sl = zeros(Nz,length(t_total));
rootC_14st = zeros(Nz,length(t_total));
%% Initial conditions bulk C 

Crz_i       = 1000;%.*exp(-7*z);%litterC_total*0.5./a0_r;%.*(D/a0_r)*D;%(litterC_total*0.6)/a0_r.*exp(-z);%5*exp(-5.*z);%1500/21.33;%(120*exp(-z/100))*dt_a/dz;%;%1.26/a0;%(24/a0)*exp(-z/(Nz-1));%(30*exp(-z/100))/(dt_a*dz);%*.001; %in mol/m^2year, concentration of C in rapid pool (veg input)
% changed from (120*exp(-z/100))*dt_a/dz
Cslz_i      = 500;%.*exp(-7*z);%0;%(6*exp(-z/100))/(dt_a*dz);%(litterC_total*0.4*dt_a/dz).*az_sl.*(D/a0_sl)*D;%1700.*exp(-z);%(litterC_total*0.4)/a0_sl.*exp(-z);%1*exp(-z); %/4; %;%27.3/a0;%(6*exp(-z/100))*dt_a/dz; %C concentration in slow pool, veg input 
Cstz_i      =  10;%.*exp(-7*z);%(litterC_total*0.1*dt_a/dz).*az_st.*(D/a0_st)*D;%Crz_i/12; %13.86/a0;%%(2*exp(-z/100))*dt_a/dz; %C concentration of steady pool

% Crz_i  = SpUp.Crz_i;
% Cslz_i = SpUp.Cslz_i;
% Cstz_i = SpUp.Cstz_i;
if SUP == 1
    Crz_i  = SpUp.Crz_i;
    Cslz_i = SpUp.Cslz_i;
    Cstz_i = SpUp.Cstz_i;
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Crz         = zeros(Nz,length(t_total));
Cslz        = zeros(Nz,length(t_total));
Cstz        = zeros(Nz,length(t_total));

Crz(:,1)    = Crz_i;
Cslz(:,1)   = Cslz_i;
Cstz(:,1)   = Cstz_i;

%% Initial Conditions 13C

Crz_13_i       = (RC(1).*Crz_i)./(1+RC(1)); %in mol/m^2year, concentration of 13C in rapid pool (veg input)
Cslz_13_i      = (RC(1).*Cslz_i)./(1+RC(1)); %13C concentration in slow pool, veg input 
Cstz_13_i      = (RC(1).*Cstz_i)./(1+RC(1)); %13C concentration of steady pool

% Crz_13_i  = SpUp.Crz_13_i;
% Cslz_13_i = SpUp.Cslz_13_i;
% Cstz_13_i = SpUp.Cstz_13_i;
if SUP == 1
    Crz_13_i  = SpUp.Crz_13_i;
    Cslz_13_i = SpUp.Cslz_13_i;
    Cstz_13_i = SpUp.Cstz_13_i;    
end 
Crz_13        = zeros(Nz,length(t_total));
Cslz_13       = zeros(Nz,length(t_total));
Cstz_13       = zeros(Nz,length(t_total));

Crz_13(:,1)    = Crz_13_i;
Cslz_13(:,1)   = Cslz_13_i;
Cstz_13(:,1)   = Cstz_13_i;

%% Intitial Conditions 12C 
Crz_12_i       = Crz_i.*(1/(RC(1)+1)); %in mol/m^2year, concentration of 12C in rapid pool (veg input)
Cslz_12_i      = Cslz_i.*(1/(RC(1)+1)); %12C concentration in slow pool, veg input 
Cstz_12_i      = Cstz_i.*(1/(RC(1)+1)); %12C concentration of steady pool

% Crz_12_i  = SpUp.Crz_12_i;
% Cslz_12_i = SpUp.Cslz_12_i;
% Cstz_12_i = SpUp.Cstz_12_i;
if SUP == 1
    Crz_12_i  = SpUp.Crz_12_i;
    Cslz_12_i = SpUp.Cslz_12_i;
    Cstz_12_i = SpUp.Cstz_12_i;    
end 
Crz_12        = zeros(Nz,length(t_total));
Cslz_12       = zeros(Nz,length(t_total));
Cstz_12       = zeros(Nz,length(t_total));

Crz_12(:,1)    = Crz_12_i;
Cslz_12(:,1)   = Cslz_12_i;
Cstz_12(:,1)   = Cstz_12_i;

%% Initial Conditions 14C 
% Same as above but done using mol convert function
Crz_14_i = Crz_i*(10e-12);%molCconvert(Aabs*dt_a,d13Crz_i,I_r_i,Crz_14_i_fmc);
Cslz_14_i = Cslz_i*(10e-12);%molCconvert(Aabs*dt_a,d13Cslz_i,I_sl_i,Cslz_14_i_fmc);
Cstz_14_i = Cstz_i*(10e-12);%molCconvert(Aabs*dt_a,d13Cstz_i,I_st_i,Cstz_14_i_fmc);

% Crz_14_i  = SpUp.Crz_14_i;
% Cslz_14_i = SpUp.Cslz_14_i;
% Cstz_14_i = SpUp.Cstz_14_i;
if SUP == 1
    Crz_14_i  = SpUp.Crz_14_i;
    Cslz_14_i = SpUp.Cslz_14_i;
    Cstz_14_i = SpUp.Cstz_14_i;    
end 
Crz_14 = zeros(Nz,length(t_total));
Cslz_14 = zeros(Nz,length(t_total));
Cstz_14 = zeros(Nz,length(t_total));

Crz_14(:,1)    = Crz_14_i;
Cslz_14(:,1)   = Cslz_14_i;
Cstz_14(:,1)   = Cstz_14_i;

DC14_r = zeros(Nz,length(t_total));
DC14_sl = zeros(Nz,length(t_total));
DC14_st = zeros(Nz,length(t_total));

DC14_r_i  = ((((Crz_14_i./Crz_i)/(Aabs)).*((0.975^2)./((1+(d13Cconvert(Crz_13_i,Crz_12_i)./1000)).^2)))-1).*1000;
DC14_sl_i = ((((Cslz_14_i./Cslz_i)/(Aabs)).*((0.975^2)./((1+(d13Cconvert(Cslz_13_i,Cslz_12_i)./1000)).^2)))-1).*1000;
DC14_st_i = ((((Cstz_14_i./Cstz_i)/(Aabs)).*((0.975^2)./((1+(d13Cconvert(Cstz_13_i,Cstz_12_i)./1000)).^2)))-1).*1000;
if SUP == 1
    DC14_r_i  = SpUp.DC14_r_i;
    DC14_sl_i = SpUp.DC14_sl_i;
    DC14_st_i = SpUp.DC14_st_i;    
end 
% DC14_r_i  = SpUp.DC14_r_i;
% DC14_sl_i = SpUp.DC14_sl_i;
% DC14_st_i = SpUp.DC14_st_i;

DC14_r(:,1) = DC14_r_i;
DC14_sl(:,1) = DC14_sl_i;
DC14_st(:,1) = DC14_st_i;

%% CO2

% CO2_conc_12(j) = (time);
% CO2_conc_13(j) = mol13_CO2_atm(time);
% CO2_conc_14(j) = mol14_CO2_atm(time); 
%                 
% CO2_conc_i = molCO2_atm(1);
% CO2_conc_12_i =mol12_CO2_atm(1);
% CO2_conc_13_i = mol13_CO2_atm(1);
% CO2_conc_14_i = mol14_CO2_atm(1);

CO2_conc     = zeros(Nz,length(t_total));
CO2_conc_13  = zeros(Nz,length(t_total));
CO2_conc_12  = zeros(Nz,length(t_total));
CO2_conc_14  = zeros(Nz,length(t_total));

% CO2_conc(:,1) = CO2_conc_i;
% CO2_conc_13(:,1) = CO2_conc_13_i;
% CO2_conc_12(:,1) = CO2_conc_12_i;
% CO2_conc_14(:,1) = CO2_conc_14_i;
% 
% CO2_lost_total = zeros(Nz,length(t_total));
% CO2_lost_total_12 = zeros(Nz,length(t_total));
% CO2_lost_total_13 = zeros(Nz,length(t_total));
% CO2_lost_total_14 = zeros(Nz,length(t_total));
% 
% d12Cdz2 = zeros(Nz,1);
% d13Cdz2 = zeros(Nz,1);
% d14Cdz2 = zeros(Nz,1);
% 
% dDs12Cdz = zeros(Nz,1);
% dDs13Cdz = zeros(Nz,1);
% dDs14Cdz = zeros(Nz,1);


%% Initializing Matrix Parameters of variables 
% length(t_total) = amount of time (endtime) over dt
% Nz is the # of nodes in the z direction (101 for 1 m of soil)

Advgradient_r = zeros(Nz,length(t_total));
Advgradient_13r = zeros(Nz,length(t_total));
Advgradient_12r = zeros(Nz,length(t_total));
Advgradient_14r = zeros(Nz,length(t_total));

Advgradient_sl = zeros(Nz,length(t_total));
Advgradient_13sl = zeros(Nz,length(t_total));
Advgradient_12sl = zeros(Nz,length(t_total));
Advgradient_14sl = zeros(Nz,length(t_total));

Advgradient_st = zeros(Nz,length(t_total));
Advgradient_13st = zeros(Nz,length(t_total));
Advgradient_12st = zeros(Nz,length(t_total));
Advgradient_14st = zeros(Nz,length(t_total));

%% Initializing diffusion
Cdiff_r = zeros(Nz,length(t_total));
Cdiff_13r = zeros(Nz,length(t_total));
Cdiff_12r = zeros(Nz,length(t_total));
Cdiff_14r = zeros(Nz,length(t_total));

Cdiff_sl = zeros(Nz,length(t_total));
Cdiff_13sl = zeros(Nz,length(t_total));
Cdiff_12sl = zeros(Nz,length(t_total));
Cdiff_14sl = zeros(Nz,length(t_total));

Cdiff_st = zeros(Nz,length(t_total));
Cdiff_13st = zeros(Nz,length(t_total));
Cdiff_12st = zeros(Nz,length(t_total));
Cdiff_14st = zeros(Nz,length(t_total));

%% Mass balance factors
% Factor to determine the proportion of 13C or 12C CO2
%these are initialized here using the d13C value of the vegetation, i.e.
%RC(1). They are then updated in the loop as the d13C value of each pool
%changes.

CO2_12fac_r = (1./(1 + (alpha_r*RC(1)))).*ones(Nz,1);
CO2_13fac_r = (alpha_r*RC(1)./(1+(alpha_r*RC(1)))).*ones(Nz,1);

CO2_12fac_sl = (1./(1 + (alpha_sl*RC(1)))).*ones(Nz,1);
CO2_13fac_sl = (alpha_sl*RC(1)./(1+(alpha_sl*RC(1)))).*ones(Nz,1);

CO2_12fac_st = (1./(1 + (alpha_st*RC(1)))).*ones(Nz,1);
CO2_13fac_st = (alpha_st*RC(1)./(1+(alpha_st*RC(1)))).*ones(Nz,1);

C14_fac = 10e-12.*ones(Nz,1);

%% Initializing C ratio 
RC_total    = zeros(Nz,length(t_total));
%RC = (Rpdb*((d13Cp/1000)+1)).*ones(1,length(t_total));
RC_i = Crz_13_i./Crz_12_i; 
RC_total(1,:) = RC; %initial condition
RC_total(:,1) = RC_i; %boundary condition

RCrz = RC_total;
RCslz = RC_total;
RCstz = RC_total;

d13Crz = ((RCrz./Rpdb) -1)*1000;
d13Cslz = d13Crz;
d13Cstz = d13Crz;

%% Initializing components in the time and depth loop 
C_consumed_r = zeros(Nz,length(t_total));
C_consumed_r_13 = zeros(Nz,length(t_total));
C_consumed_r_12 = zeros(Nz,length(t_total));
C_consumed_r_14 = zeros(Nz,length(t_total));

C_consumed_sl = zeros(Nz,length(t_total));
C_consumed_sl_13 = zeros(Nz,length(t_total));
C_consumed_sl_12 = zeros(Nz,length(t_total));
C_consumed_sl_14 = zeros(Nz,length(t_total));

C_consumed_st = zeros(Nz,length(t_total));
C_consumed_st_13 = zeros(Nz,length(t_total));
C_consumed_st_12 = zeros(Nz,length(t_total));
C_consumed_st_14 = zeros(Nz,length(t_total));
                
CO2_lost_r   = zeros(Nz,length(t_total));
CO2_lost_13r = zeros(Nz,length(t_total));
CO2_lost_12r = zeros(Nz,length(t_total));
CO2_lost_14r = zeros(Nz,length(t_total));

CO2_lost_sl   = zeros(Nz,length(t_total));
CO2_lost_13sl = zeros(Nz,length(t_total));
CO2_lost_12sl = zeros(Nz,length(t_total));
CO2_lost_14sl = zeros(Nz,length(t_total));

CO2_lost_st  = zeros(Nz,length(t_total));
CO2_lost_13st = zeros(Nz,length(t_total));
CO2_lost_12st = zeros(Nz,length(t_total));
CO2_lost_14st  = zeros(Nz,length(t_total));

growth_r = zeros(Nz,length(t_total));
transformed_C_r_to_sl = zeros(Nz,length(t_total));

growth_r_13 = zeros(Nz,length(t_total));
growth_r_12 = zeros(Nz,length(t_total));
growth_r_14 = zeros(Nz,length(t_total)); 

transformed_C_r_to_sl_13 = zeros(Nz,length(t_total));
transformed_C_r_to_sl_12 = zeros(Nz,length(t_total));
transformed_C_r_to_sl_14 = zeros(Nz,length(t_total));

growth_sl   =  zeros(Nz,length(t_total));
transformed_C_sl_to_st = zeros(Nz,length(t_total));

growth_sl_13 = zeros(Nz,length(t_total));
growth_sl_12 = zeros(Nz,length(t_total));
growth_sl_14 = zeros(Nz,length(t_total));

transformed_C_sl_to_st_13 = zeros(Nz,length(t_total));
transformed_C_sl_to_st_12 = zeros(Nz,length(t_total));
transformed_C_sl_to_st_14 = zeros(Nz,length(t_total));
% 
CO2_lost_total = zeros(Nz,length(t_total)); 
CO2_lost_total_12 = zeros(Nz,length(t_total));
CO2_lost_total_13 = zeros(Nz,length(t_total));
CO2_lost_total_14 = zeros(Nz,length(t_total));