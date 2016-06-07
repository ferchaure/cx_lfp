function par = par_cb_lfp() 
 
par.IP_address = '': %string, empty for default mode 
par.channels = 97:115;
 
%periodogram 
par.chunk_time = 1.1; %seconds
par.fmax_update = 3000; 
 
%plotting 
par.fmin_disp = 1; 
par.fmax_disp= 3000; 
 
%power fit 
par.fini_fit = 2; 
par.fmid_fit = 300; 
end 
