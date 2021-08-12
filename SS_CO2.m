function [CO2_old,CO2_old_12,CO2_old_13,CO2_old_14] = SS_CO2(time,dz,Nz,D_diff,poros,dt_diff, Aabs, d13Ca, pCO2, C14_atm, Ratm,CO2_lost_total,CO2_lost_total_12,CO2_lost_total_13,CO2_lost_total_14,...
                      CO2_conc,CO2_conc_12,CO2_conc_13,CO2_conc_14)
    % Atmospheric values in moles              
    molCO2_atm(time)    = pCO2(time)/.022400/10^6;
    mol12_CO2_atm(time) = molCO2_atm(time)/(Ratm(time)+1);
    mol13_CO2_atm(time) = (Ratm(time)*molCO2_atm(time))/(Ratm(time)+1);
    mol14_CO2_atm(time) = (((C14_atm(time)/1000)+1)/((0.975^2)/((1+(d13Ca(time)/1000))^2)))*Aabs*molCO2_atm(time);
    
    D_12C= D_diff;
    D_13C =D_diff/1.0044;
    D_14C =D_diff/1.00868; % W Yang 1994
    fart_time = 0.005;
    
%%
   if time == 1 % When time = 1, the initial CO2 profile is that of atmospheric values 
            for p = 1:Nz
                CO2_conc(p)    = molCO2_atm(time);
                CO2_conc_12(p) = mol12_CO2_atm(time);
                CO2_conc_13(p) = mol13_CO2_atm(time);
                CO2_conc_14(p) = mol14_CO2_atm(time); 
            end 
            CO2_old(1:Nz) = CO2_conc(1:Nz);
            CO2_old_12(1:Nz) = CO2_conc_12(1:Nz);
            CO2_old_13(1:Nz) = CO2_conc_13(1:Nz);
            CO2_old_14(1:Nz) = CO2_conc_14(1:Nz);
            
     CO2_diff = 1;     
     while CO2_diff > 0.1   
            
            for p = 2:Nz-1
                dCdz2 (p)   = (CO2_old(p-1) + CO2_old(p+1) - 2*CO2_old(p) )/(dz^2); %mol/m5
                d12Cdz2 (p) = (CO2_old_12(p-1) + CO2_old_12(p+1) - 2*CO2_old_12(p) )/(dz^2);
                d13Cdz2(p)  = (CO2_old_13(p-1) + CO2_old_13(p+1) - 2*CO2_old_13(p) )/(dz^2);
                d14Cdz2(p)  = (CO2_old_14(p-1) + CO2_old_14(p+1) - 2*CO2_old_14(p) )/(dz^2);
                
                delC(p)   = dt_diff*((D_diff(p)*dCdz2(p))+(CO2_lost_total(p,time)/poros));
                del12C(p) = dt_diff*((D_12C(p)*d12Cdz2 (p))+(CO2_lost_total_12(p,time)/poros));
                del13C(p) = dt_diff*((D_13C(p)*d13Cdz2(p))+(CO2_lost_total_13(p,time)/poros));
                del14C(p) = dt_diff*((D_14C(p)*d14Cdz2(p))+(CO2_lost_total_14(p,time)/poros));
                
                CO2_new(p)    = delC(p) + CO2_old(p);
                CO2_new_12(p) = del12C(p) + CO2_old_12(p);
                CO2_new_13(p) = del13C(p) + CO2_old_13(p);
                CO2_new_14(p) = del14C(p) + CO2_old_14(p);

            end 
                CO2_new(1) = molCO2_atm(time);
                CO2_new_12(1) = mol12_CO2_atm(time);
                CO2_new_13(1) =  mol13_CO2_atm(time);
                CO2_new_14(1) = mol14_CO2_atm(time);

                CO2_new(Nz) = CO2_new(Nz-1);
                CO2_new_12(Nz) = CO2_new_12(Nz-1);
                CO2_new_13(Nz) =  CO2_new_13(Nz-1);
                CO2_new_14(Nz) = CO2_new_14(Nz-1);
            for p = 2:Nz-1
                dCdz2 (p) = (CO2_new(p-1) + CO2_new(p+1) - 2*CO2_new(p) )/(dz^2); %mol/m5
                d12Cdz2 (p) = (CO2_new_12(p-1) + CO2_new_12(p+1) - 2*CO2_new_12(p) )/(dz^2);
                d13Cdz2(p) = (CO2_new_13(p-1) + CO2_new_13(p+1) - 2*CO2_new_13(p) )/(dz^2);
                d14Cdz2(p) = (CO2_new_14(p-1) + CO2_new_14(p+1) - 2*CO2_new_14(p) )/(dz^2);
                
                delC_new(p)   = dt_diff*((D_diff(p)*dCdz2(p))+(CO2_lost_total(p,time)/poros));
                del12C_new(p) = dt_diff*((D_12C(p)*d12Cdz2 (p))+(CO2_lost_total_12(p,time)/poros));
                del13C_new(p) = dt_diff*((D_13C(p)*d13Cdz2(p))+(CO2_lost_total_13(p,time)/poros));
                del14C_new(p) = dt_diff*((D_14C(p)*d14Cdz2(p))+(CO2_lost_total_14(p,time)/poros));
                
                CO2_new(p)    = CO2_old(p) + (delC(p) + delC_new(p))/2;
                CO2_new_12(p) = CO2_old_12(p) + (del12C(p) + del12C_new(p))/2;
                CO2_new_13(p) = CO2_old_13(p) + (del13C(p) + del13C_new(p))/2;
                CO2_new_14(p) = CO2_old_14(p) + (del14C(p) + del14C_new(p))/2;
            end 
                CO2_new(1) = molCO2_atm(time);
                CO2_new_12(1) = mol12_CO2_atm(time);
                CO2_new_13(1) =  mol13_CO2_atm(time);
                CO2_new_14(1) = mol14_CO2_atm(time);

                CO2_new(Nz) = CO2_new(Nz-1);
                CO2_new_12(Nz) = CO2_new_12(Nz-1);
                CO2_new_13(Nz) =  CO2_new_13(Nz-1);
                CO2_new_14(Nz) = CO2_new_14(Nz-1);
    
            if norm(CO2_new(1:Nz) - CO2_old(1:Nz))/norm(CO2_new(1:Nz)) < fart_time
                CO2_diff = 0;
            end 
%             CO2_diff = abs(mean(CO2_new(1:Nz) - CO2_old(1:Nz)));
     
            CO2_old(1:Nz) = CO2_new(1:Nz);
            CO2_old_12(1:Nz) = CO2_new_12(1:Nz);
            CO2_old_13(1:Nz) = CO2_new_13(1:Nz);
            CO2_old_14(1:Nz) = CO2_new_14(1:Nz);
%             figure(1)
%             plot(CO2_old(1:Nz),(-1*(1:Nz)))
%             hold on
     end   
     
        CO2_conc(:,time) = CO2_old(1:Nz);
        CO2_conc_12(:,time) = CO2_old_12(1:Nz);
        CO2_conc_13(:,time) = CO2_old_13(1:Nz);
        CO2_conc_14(:,time) = CO2_old_14(1:Nz);
       
%    end        
   elseif time > 1
            
 %%
%             for j = 1:Nz
%                 CO2_conc(j)    = CO2_conc(j,time);
%                 CO2_conc_12(j) = CO2_conc_12(j,time);
%                 CO2_conc_13(j) = CO2_conc_13(j,time);
%                 CO2_conc_14(j) = CO2_conc_14(j,time); 
%             end 
            CO2_old(1:Nz) = CO2_conc(1:Nz);
            CO2_old_12(1:Nz) = CO2_conc_12(1:Nz);
            CO2_old_13(1:Nz) = CO2_conc_13(1:Nz);
            CO2_old_14(1:Nz) = CO2_conc_14(1:Nz);
            
     CO2_diff = 1;     
     while CO2_diff > 0.1    
            
            for p = 2:Nz-1
                dCdz2 (p)   = (CO2_old(p-1) + CO2_old(p+1) - 2*CO2_old(p) )/(dz^2); %mol/m5
                d12Cdz2 (p) = (CO2_old_12(p-1) + CO2_old_12(p+1) - 2*CO2_old_12(p) )/(dz^2);
                d13Cdz2(p)  = (CO2_old_13(p-1) + CO2_old_13(p+1) - 2*CO2_old_13(p) )/(dz^2);
                d14Cdz2(p)  = (CO2_old_14(p-1) + CO2_old_14(p+1) - 2*CO2_old_14(p) )/(dz^2);
                
                delC(p)   = dt_diff*((D_diff(p)*dCdz2(p))+(CO2_lost_total(p,time)/poros));
                del12C(p) = dt_diff*((D_12C(p)*d12Cdz2 (p))+(CO2_lost_total_12(p,time)/poros));
                del13C(p) = dt_diff*((D_13C(p)*d13Cdz2(p))+(CO2_lost_total_13(p,time)/poros));
                del14C(p) = dt_diff*((D_14C(p)*d14Cdz2(p))+(CO2_lost_total_14(p,time)/poros));
                
                CO2_new(p)    = delC(p) + CO2_old(p);
                CO2_new_12(p) = del12C(p) + CO2_old_12(p);
                CO2_new_13(p) = del13C(p) + CO2_old_13(p);
                CO2_new_14(p) = del14C(p) + CO2_old_14(p);
                
            end 
                CO2_new(1) = molCO2_atm(time);
                CO2_new_12(1) = mol12_CO2_atm(time);
                CO2_new_13(1) =  mol13_CO2_atm(time);
                CO2_new_14(1) = mol14_CO2_atm(time);

                CO2_new(Nz) = CO2_new(Nz-1);
                CO2_new_12(Nz) = CO2_new_12(Nz-1);
                CO2_new_13(Nz) =  CO2_new_13(Nz-1);
                CO2_new_14(Nz) = CO2_new_14(Nz-1);
                
            for p = 2:Nz-1
                dCdz2 (p) = (CO2_new(p-1) + CO2_new(p+1) - 2*CO2_new(p) )/(dz^2); %mol/m5
                d12Cdz2 (p) = (CO2_new_12(p-1) + CO2_new_12(p+1) - 2*CO2_new_12(p) )/(dz^2);
                d13Cdz2(p) = (CO2_new_13(p-1) + CO2_new_13(p+1) - 2*CO2_new_13(p) )/(dz^2);
                d14Cdz2(p) = (CO2_new_14(p-1) + CO2_new_14(p+1) - 2*CO2_new_14(p) )/(dz^2);

                delC_new(p)   = dt_diff*((D_diff(p)*dCdz2(p))+(CO2_lost_total(p,time)/poros));
                del12C_new(p) = dt_diff*((D_12C(p)*d12Cdz2 (p))+(CO2_lost_total_12(p,time)/poros));
                del13C_new(p) = dt_diff*((D_13C(p)*d13Cdz2(p))+(CO2_lost_total_13(p,time)/poros));
                del14C_new(p) = dt_diff*((D_14C(p)*d14Cdz2(p))+(CO2_lost_total_14(p,time)/poros));
                
                CO2_new(p)    = CO2_old(p) + (delC(p) + delC_new(p))/2;
                CO2_new_12(p) = CO2_old_12(p) + (del12C(p) + del12C_new(p))/2;
                CO2_new_13(p) = CO2_old_13(p) + (del13C(p) + del13C_new(p))/2;
                CO2_new_14(p) = CO2_old_14(p) + (del14C(p) + del14C_new(p))/2;
            end 
                CO2_new(1) = molCO2_atm(time);
                CO2_new_12(1) = mol12_CO2_atm(time);
                CO2_new_13(1) =  mol13_CO2_atm(time);
                CO2_new_14(1) = mol14_CO2_atm(time);

                CO2_new(Nz) = CO2_new(Nz-1);
                CO2_new_12(Nz) = CO2_new_12(Nz-1);
                CO2_new_13(Nz) =  CO2_new_13(Nz-1);
                CO2_new_14(Nz) = CO2_new_14(Nz-1);
                
            
            if norm(CO2_new(1:Nz) - CO2_old(1:Nz))/norm(CO2_new(1:Nz)) < fart_time
                CO2_diff = 0;
            end 
     
            CO2_old(1:Nz) = CO2_new(1:Nz);
            CO2_old_12(1:Nz) = CO2_new_12(1:Nz);
            CO2_old_13(1:Nz) = CO2_new_13(1:Nz);
            CO2_old_14(1:Nz) = CO2_new_14(1:Nz); 
            
%             figure(1)
%             plot(CO2_new(1:Nz))
%             hold on
     end   
     
        CO2_conc(:,time) = CO2_old(1:Nz);
        CO2_conc_12(:,time) = CO2_old_12(1:Nz);
        CO2_conc_13(:,time) = CO2_old_13(1:Nz);
        CO2_conc_14(:,time) = CO2_old_14(1:Nz);
            
            
            
        end    
end 