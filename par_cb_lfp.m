function par = par_cb_lfp() 
%just for one nsp
par.IP_address = ''; %string, empty for default mode 
par.channels = 'all'; %vector of channels ID or all to use all channels
par.ffftlength = 2^ceil(log2(30000*2))/30000; %seconds it will be increased to use a power of 2 
par.freq_line = 60; 
par.notch_width = 1;

%filters
par.num_notchs = 50;

par.custom_filter = struct;
par.custom_filter.('f500').enable = false;
par.custom_filter.('f500').bp1 = 2;
par.custom_filter.('f500').bp2 =  100;
par.custom_filter.('f1000').enable = false;
par.custom_filter.('f1000').bp1 = 2;
par.custom_filter.('f1000').bp2 =  200;
par.custom_filter.('f2000').enable = false;
par.custom_filter.('f2000').bp1 = 2;
par.custom_filter.('f2000').bp2 =  500;
par.custom_filter.('f10000').enable = false;
par.custom_filter.('f10000').bp1 = 2;
par.custom_filter.('f10000').bp2 =  3000;
par.custom_filter.('f30000').enable = false;
par.custom_filter.('f30000').bp1 = 2;
par.custom_filter.('f30000').bp2 =  3000;

%plotting
par.n_blocks  = 2;

par.x_power_manual = struct;
par.x_power_manual.('f500').min = 0;
par.x_power_manual.('f500').max =  250;
par.x_power_manual.('f500').min_zoom = 0;
par.x_power_manual.('f500').max_zoom = 200;

par.x_power_manual.('f1000').min = 0;
par.x_power_manual.('f1000').max =  500;
par.x_power_manual.('f1000').min_zoom = 0;
par.x_power_manual.('f1000').max_zoom = 200;

par.x_power_manual.('f2000').min = 0;
par.x_power_manual.('f2000').max =  1000;
par.x_power_manual.('f2000').min_zoom = 0;
par.x_power_manual.('f2000').max_zoom = 200;

par.x_power_manual.('f10000').min = 0;
par.x_power_manual.('f10000').max =  3000;
par.x_power_manual.('f10000').min_zoom = 0;
par.x_power_manual.('f10000').max_zoom = 300;

par.x_power_manual.('f30000').min = 0;
par.x_power_manual.('f30000').max =  3000;
par.x_power_manual.('f30000').min_zoom = 0;
par.x_power_manual.('f30000').max_zoom = 300;

%extra

par.fstop_h = 1.3;%times fmax of bandpass (bp2) this have to be calculated in a gui
par.fstop_l = 0.5; 
par.Rp = 20;
par.Rs = 0.07;
end
 
