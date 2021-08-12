        
function interp14C = C14fit(ntotal,run_time,end_date)

%ntotal: total number of time iterations
%y1: endtime (e.g. 500 years)
%y2: start time (0)

%C14data = xlsread('radiocarbon_data_all.xlsx'); %from -61 to 50000 years BP
% load C14data3
load C14data5
% [~,start_time] = min(abs(y1-C14data2(:,4)));
% [~,end_time]   = min(abs(y2-C14data2(:,4))); %3rd column is calendar year and ending in 2011
[~,end_time]   = min(abs(end_date-C14data5(:,3))); % Where the last data point is/where to stop interpolating 
[~,start_time] = min(abs((end_date-run_time)-C14data5(:,3))); 

newtimes_14C = linspace(C14data5(end_time,4),C14data5(start_time,4),ntotal);
C14interp = interp1(C14data5(end_time:start_time,4),...
                C14data5(end_time:start_time,2),newtimes_14C);

interp14C = C14interp;

end 

    