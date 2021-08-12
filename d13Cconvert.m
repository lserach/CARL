function d13C = d13Cconvert(C13,C12)
Rpdb = 0.0112372;
d13C = (((C13./C12)./Rpdb)-1).*1000;
end 