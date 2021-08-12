        
function interpd13C = d13Cfit(ntotal,run_time,end_date)
 
%ntotal: total number of time iterations
%y1: the start time (yr)
%y2: the end time (yr)
% load d13Cdata2
load d13Cdata3
%d13Cdata = xlsread('d13C_data_all.xlsx');

% [~,end_time]   = min(abs(y2-d13Cdata(:,4))); % Where the last data point is/where to stop interpolating 
% [~,start_time] = min(abs(endtime-d13Cdata(:,4))); 

[~,end_time]   = min(abs(end_date-d13Cdata3(:,3))); % Where the last data point is/where to stop interpolating 
[~,start_time] = min(abs((end_date-run_time)-d13Cdata3(:,3))); 


newtimes_d13C = linspace(d13Cdata3((end_time),4),d13Cdata3((start_time + 1),4),ntotal);%linspace(y1,y2,ntotal);
d13Cinterp = interp1(d13Cdata3(end_time:(start_time + 1),4),...
                     d13Cdata3(end_time:(start_time + 1),2),...
                     newtimes_d13C);

interpd13C = d13Cinterp;



end 