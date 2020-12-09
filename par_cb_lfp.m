function par = par_cb_lfp() 
%just for one nsp
par.IP_address = '192.168.137.3'; %string, empty for default mode 
par.channels = 'all'; %vector of channels ID or all to use all channels
par.ffftlength = 2^ceil(log2(30000*2))/30000; %seconds it will be increased to use a power of 2 
%^this has to be the largest one in the used freqs
par.freqs_notch = [60 120]; 
par.notch_width = 1;
%plotting
par.n_blocks  = 5;

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

%filters
par.custom_filter = struct;
par.custom_filter.('f500').enable = false;
par.custom_filter.('f500').bp1 = 2;
par.custom_filter.('f500').bp2 =  100;
par.custom_filter.('f500').order =  2;
par.custom_filter.('f500').filter_type =  'ellip_order'; %'ellip_order' or 'ellip_stop_band'
par.custom_filter.('f500').fstop1 = 0.5;
par.custom_filter.('f500').fstop2 = 130;
par.custom_filter.('f500').Rp = 0.07;
par.custom_filter.('f500').Rs = 20;

par.custom_filter.('f1000').enable = false;
par.custom_filter.('f1000').bp1 = 2;
par.custom_filter.('f1000').bp2 =  200;
par.custom_filter.('f1000').enable = false;
par.custom_filter.('f1000').filter_type =  'ellip_order'; %'ellip_order' or 'ellip_stop_band'
par.custom_filter.('f1000').fstop1 = 0.5;
par.custom_filter.('f1000').fstop2 = 260;
par.custom_filter.('f1000').Rp = 0.07;
par.custom_filter.('f1000').Rs = 20;

par.custom_filter.('f2000').enable = false;
par.custom_filter.('f2000').bp1 = 2;
par.custom_filter.('f2000').bp2 =  500;
par.custom_filter.('f2000').order =  2;
par.custom_filter.('f2000').filter_type =  'ellip_order'; %'ellip_order' or 'ellip_stop_band'
par.custom_filter.('f2000').fstop1 = 0.5;
par.custom_filter.('f2000').fstop2 = 800;
par.custom_filter.('f2000').Rp = 0.07;
par.custom_filter.('f2000').Rs = 20;

par.custom_filter.('f10000').enable = false;
par.custom_filter.('f10000').bp1 = 2;
par.custom_filter.('f10000').bp2 =  3000;
par.custom_filter.('f10000').order =  2;
par.custom_filter.('f10000').filter_type =  'ellip_order'; %'ellip_order' or 'ellip_stop_band'
par.custom_filter.('f10000').fstop1 = 0.5;
par.custom_filter.('f10000').fstop2 = 4000;
par.custom_filter.('f10000').Rp = 0.07;
par.custom_filter.('f10000').Rs = 20;

par.custom_filter.('f30000').enable = 1;
par.custom_filter.('f30000').bp1 = 300;
par.custom_filter.('f30000').bp2 =  3000;
par.custom_filter.('f30000').order =  2;
par.custom_filter.('f30000').filter_type =  'ellip_order'; %'ellip_order' or 'ellip_stop_band'
par.custom_filter.('f30000').fstop1 = 0.5;
par.custom_filter.('f30000').fstop2 = 4000;
par.custom_filter.('f30000').Rp = 0.07;
par.custom_filter.('f30000').Rs = 20;

end
 
