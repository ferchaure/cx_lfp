function par = par_cb_lfp() 
 
par.channels = 97:115; 
 
%periodogram 
par.n_samples = 2^15; 
par.fmax_update = 3000; 
 
%plotting 
par.fmin_disp = 1; 
par.fmax_disp= 3000; 
 
%power fit 
par.fini_fit = 2; 
par.fmid_fit = 300; 
end 
