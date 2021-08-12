function interpppmV = pCO2fit(ntotal,run_time,end_date)

%ntotal: total number of time iterations
%y1: endtime (e.g. 500 years)
%y2: start time (0)

%ppmVdata = xlsread('pCO2_data.xlsx');
% load ppmVdata
load ppmVdata2

% [~,start_time] = min(abs(y1-ppmVdata(:,4))); %min(abs(y1-ppmVdata(:,3))); 
% [~,end_time]   = min(abs(y2-ppmVdata(:,4))); %3rd column is calendar year and ending in 2011
[~,end_time]   = min(abs(end_date-ppmVdata2(:,3))); % Where the last data point is/where to stop interpolating 
[~,start_time] = min(abs((end_date-run_time)-ppmVdata2(:,3))); 

% 4th column is years BP (present is 2011) min(abs(y2-ppmVdata(:,3)));

newtimes_ppmV = linspace(ppmVdata2(end_time,4),ppmVdata2(start_time,4),ntotal); %linspace(y1,y2,ntotal);

ppmVinterp = interp1(ppmVdata2(end_time:(start_time +1) ,4),...
                     ppmVdata2(end_time:(start_time +1) ,2),...
                     newtimes_ppmV);

interpppmV = ppmVinterp;

end
 
    