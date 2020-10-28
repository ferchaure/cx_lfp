function varargout = cx_lfp(varargin)
% cx_lfp MATLAB code for cx_lfp.fig
%      cx_lfp, by itself, creates a new cx_lfp or raises the existing
%      singleton*.
%
%      H = cx_lfp returns the handle to a new cx_lfp or the handle to
%      the existing singleton*.
%
%      cx_lfp('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in cx_lfp.M with the given input arguments.
%
%      cx_lfp('Property','Value',...) creates a new cx_lfp or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cx_lfp_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cx_lfp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cx_lfp

% Last Modified by GUIDE v2.5 27-Oct-2020 17:19:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cx_lfp_OpeningFcn, ...
                   'gui_OutputFcn',  @cx_lfp_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
end

% --- Executes just before cx_lfp is made visible.
function cx_lfp_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cx_lfp (see VARARGIN)

% Choose default command line output for cx_lfp
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
start_adq(hObject)
end

% --- Outputs from this function are returned to the command line.
function varargout = cx_lfp_OutputFcn(~, ~, handles) 
varargout{1} = handles.output;
fig = gcf;
fig.WindowState='maximized';
end

function display_channel(n,hObject)
    h = guidata(hObject);
    set(h.N_lb,'String',num2str(h.N));
    set(h.channels_lb,'Value',n);

    if strcmp(h.ch_data{n}.unit,'uV') && h.par.custom_filter.(h.ch_data{n}.parsr).enable
        additional_curve = 'Notch and Custom Filter'; 
    else
        additional_curve = 'Notch'; 
    end
    log_psd = 10*log10(h.ch_data{n}.psd);
    log_psd_filtered = 10*log10(h.ch_data{n}.psd_filtered);
    set(h.channel_label,'String',sprintf('%s | ID: %d', h.ch_data{n}.label,h.ch_data{n}.ch))

    %set gui params
    set(h.info_channel,'String',{sprintf('Sampl. Rate: %1.0f Hz',h.ch_data{n}.sr),...
        sprintf('Cont. Acq.: %s',h.ch_data{n}.smpfilter)});
    
    pars = h.ch_data{n}.parsr;
    set(h.xpower_min,'String',num2str(h.par.x_power_manual.(pars).min));
    set(h.xpower_max,'String',num2str(h.par.x_power_manual.(pars).max));
    set(h.xpower_z_min,'String',num2str(h.par.x_power_manual.(pars).min_zoom));
    set(h.xpower_z_max,'String',num2str(h.par.x_power_manual.(pars).max_zoom));  
    
    gui_par = [h.ch_data{n}.parsr  h.ch_data{n}.unit];
    
    ylim_el = {'ypower_min','ypower_max','ypower_z_min','ypower_z_max','yraw_min','yraw_max','yfiltered_min','yfiltered_max'};
    ylim_cv = {'fix_ypower_z_cb','fix_ypower_cb','fix_yraw_cb','fix_yfiltered_cb'};
    for ed  = ylim_el
        set(h.(ed{1}),'String',num2str(getappdata(h.(ed{1}),gui_par)))
    end
    for cb  = ylim_cv
        v = getappdata(h.(cb{1}),gui_par);
        set(h.(cb{1}),'Value',v);
        enable_editlines(h.(cb{1}),v,h)
    end
    %spectrum
    cla(h.spectrum); hold(h.spectrum,'on');
    plot(h.spectrum,h.sr_data{h.ch_data{n}.sri}.fs,log_psd,'LineWidth',1.5);
    if getappdata(h.fix_ypower_cb,gui_par)
        yl = [getappdata(h.ypower_min,gui_par),getappdata(h.ypower_max,gui_par)];
    else
        ylim(h.spectrum,[-inf, inf]);
        yl = ylim(h.spectrum);
    end
    plot(h.spectrum,h.sr_data{h.ch_data{n}.sri}.fs,log_psd_filtered,'r','LineWidth',1.5)
    ylim(h.spectrum,yl)
    h.spectrum.XMinorGrid = 'on'; h.spectrum.YMinorGrid = 'on';    
    xlim(h.spectrum,[h.par.x_power_manual.(pars).min,h.par.x_power_manual.(pars).max]);
    xlabel(h.spectrum,'Frequency (Hz)')
    ylabel(h.spectrum,'Power Spectrum (dB/Hz)')
    legend(h.spectrum,'Raw',additional_curve)
    
    %spectrum_zoom
    cla(h.spectrum_zoom); hold(h.spectrum_zoom,'on');
    plot(h.spectrum_zoom,h.sr_data{h.ch_data{n}.sri}.fs,log_psd,'LineWidth',1.5);
    if getappdata(h.fix_ypower_z_cb,gui_par)
        yl = [getappdata(h.ypower_z_min,gui_par),getappdata(h.ypower_z_max,gui_par)];
    else
        ylim(h.spectrum_zoom,[-inf, inf]);
        yl = ylim(h.spectrum_zoom);
    end
    plot(h.spectrum_zoom,h.sr_data{h.ch_data{n}.sri}.fs,log_psd_filtered,'r','LineWidth',1.5)
    ylim(h.spectrum_zoom,yl)
    h.spectrum_zoom.XMinorGrid = 'on'; h.spectrum_zoom.YMinorGrid = 'on';
    xlim(h.spectrum_zoom,[h.par.x_power_manual.(pars).min_zoom,h.par.x_power_manual.(pars).max_zoom]);
    xlabel(h.spectrum_zoom,'Frequency (Hz)')
    ylabel(h.spectrum_zoom,'Power Spectrum (dB/Hz)')

    %time_raw
    cla(h.time_raw)
    plot(h.time_raw,h.sr_data{h.ch_data{n}.sri}.t,h.ch_data{n}.cont,'LineWidth',1.2);
    h.time_raw.XMinorGrid = 'on'; h.time_raw.YMinorGrid = 'on';
    ylabel(h.time_raw,['Raw (' h.ch_data{n}.unit ')'])
    xlabel(h.time_raw,'Time (sec)')
    xlim(h.time_raw,[h.sr_data{h.ch_data{n}.sri}.t(1) h.sr_data{h.ch_data{n}.sri}.t(end)])
    if getappdata(h.fix_yraw_cb,gui_par)
        ylim(h.time_raw,[getappdata(h.yraw_min,gui_par),getappdata(h.yraw_max,gui_par)])
    end
    %time_filtered
    cla(h.time_filtered)
    plot(h.time_filtered,h.sr_data{h.ch_data{n}.sri}.t,h.ch_data{n}.cont_filtered,'r','LineWidth',1.2);
    h.time_raw.XMinorGrid = 'on'; h.time_filtered.YMinorGrid = 'on';
    ylabel(h.time_filtered,{additional_curve,['(' h.ch_data{n}.unit ')']})
    xlabel(h.time_filtered,'Time (sec)')
    
    
    if getappdata(h.fix_yfiltered_cb,gui_par)
        ylim(h.time_filtered,[getappdata(h.yfiltered_min,gui_par),getappdata(h.yfiltered_max,gui_par)])
    end
    xlim(h.time_filtered,[h.sr_data{h.ch_data{n}.sri}.t(1) h.sr_data{h.ch_data{n}.sri}.t(end)])
    linkaxes([h.time_raw,h.time_filtered],'x')
end


function calc_loop(~, ~, hObject)
    h = guidata(hObject);
    if h.new_data_loaded
        par = h.par;
        h.N = h.N + par.n_blocks; 
        N = h.N;
        for i = 1:length(par.channels)
            h.ch_data{i}.cont(:) =  double(h.ch_data{i}.buffer)*h.ch_data{i}.conversion;
            h.ch_data{i}.nupdated = 0;
            si = h.ch_data{i}.sri ;  %sample rate
            [psd ,~] = pwelch(h.ch_data{i}.cont,h.sr_data{si}.win,0,[],h.sr_data{si}.sr,'onesided');
            h.ch_data{i}.psd(:) = psd*par.n_blocks/N + h.ch_data{i}.psd*(N-1)/N;

            if strcmp(h.ch_data{i}.unit,'uV') && par.custom_filter.(h.ch_data{i}.parsr).enable  
                s = h.sr_data{si}.custom_filter.S;
                g = h.sr_data{si}.custom_filter.G;
            else
                s = h.sr_data{si}.notch.S;
                g = h.sr_data{si}.notch.G;
            end
            h.ch_data{i}.cont_filtered(:)= filtfilt(s,g,h.ch_data{i}.cont);
            
            [psd ,~] = pwelch(h.ch_data{i}.cont_filtered,h.sr_data{si}.win,0,[],h.sr_data{si}.sr,'onesided');
            h.ch_data{i}.psd_filtered(:) = psd*par.n_blocks/N + h.ch_data{i}.psd_filtered*(N-1)/N;       
               
        end
        h.new_data_loaded = false;
        guidata(hObject,h);
        if h.stop_refresh_cb.Value==false
            display_channel(h.channels_lb.Value,hObject)
        end
    end
end


function buffer_loop(~, ~,hObject)
%function that is called using a timer and load new data
    handles = guidata(hObject);
    %max datapoints per channel 102400
    if handles.new_data_loaded == false
        [~, x] = cbmex('trialdata',1);
        keep_loading = false;
        for ci = 1:length(handles.par.channels)
            news = size(x{ci,3},1);
            if news == 102400
                warning('lossing continuos data')
            end
            nbuffer = handles.ch_data{ci}.nupdated;
            to_update = min(news,handles.ch_data{ci}.nmax-nbuffer);
            handles.ch_data{ci}.buffer(nbuffer+1:nbuffer+to_update) = x{ci,3}(1:to_update);
            handles.ch_data{ci}.nupdated = nbuffer + to_update;
            if (handles.ch_data{ci}.nmax - handles.ch_data{ci}.nupdated)>0
                keep_loading = true;
            end
        end
        if keep_loading == false %all buffers are full
            handles.new_data_loaded = true;
        end
        guidata(hObject,handles);
    else
        [~, ~] = cbmex('trialdata',1);
    end
end
    


% --- Executes on button press in restart_button.
function restart_button_Callback(hObject, ~, ~)
% hObject    handle to restart_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)  
    stop_adq(hObject)
    start_adq(hObject)
end

% --- Executes on button press in set_param_button.
function set_param_button_Callback(hObject, eventdata, handles)
% hObject    handle to set_param_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
restart = set_par_ui();
if restart
    stop_adq(hObject)
    start_adq(hObject)
end
end

% --- Executes when user attempts to close cx_lfp.
function cx_lfp_CloseRequestFcn(hObject, eventdata, handles)
    stop_adq(hObject)
    delete(hObject);
end

function stop_adq(hObject)
    handles = guidata(hObject);
    try
        stop(handles.timer_buffer_loop);
        stop(handles.timer_calc_loop);
        cbmex('close');
        delete(handles.timer_calc_loop);
        delete(handles.timer_buffer_loop);
    catch
        warning('closed with errors')
    end
    guidata(hObject, handles);
end
    
    
function start_adq(hObject)
    if exist('C:\Program Files\Blackrock Microsystems\NeuroPort Windows Suite', 'dir')
        addpath('C:\Program Files\Blackrock Microsystems\NeuroPort Windows Suite')
    elseif exist('C:\Program Files (x86)\Blackrock Microsystems\Cerebus Windows Suite', 'dir')
        addpath('C:\Program Files (x86)\Blackrock Microsystems\Cerebus Windows Suite')
    elseif exist('C:\Program Files (x86)\Blackrock Microsystems\NeuroPort Windows Suite', 'dir')
        addpath('C:\Program Files (x86)\Blackrock Microsystems\NeuroPort Windows Suite')
    else
        warning('Using cbmex simulator')
        addpath([pwd filesep 'sim'])
    end
    
    clear functions % reset functions, force to reload set_parameters next
    handles = guidata(hObject);
    handles.par = par_cb_lfp();
    if isempty(handles.par.IP_address)
		[connection source] = cbmex('open');
	else
		[connection source] = cbmex('open', 2, 'central-addr',handles.par.IP_address, 'instance', 1);
	end
    
    handles.new_data_loaded = false;
    handles.N = 0;
    handles.n_show = 1 ;
    chs = handles.par.channels;
    if ~isstr(chs)
        cbmex('mask',0,0)
        for j = chs
            cbmex('mask',j,1)
        end
    else
        cbmex('mask',0,1)
    end
    
    cbmex('trialconfig',1,'noevent');
    [t_ini_buffer x] = cbmex('trialdata',1);
    while size(x,1) == 0
        [t_ini_buffer x] = cbmex('trialdata',1);
    end
    chs = cell2mat(x(:,1));
    handles.par.channels = chs;
    handles.cbmex_lfp.PaperPosition = handles.cbmex_lfp.Position;
    handles.cbmex_lfp.PaperPositionMode = 'auto';
    
    labels = cbmex('chanlabel',chs);
    set(handles.channels_lb,'string',labels(:,1)); 
    aux_info = cbmex('config',chs);

    handles.ch_data = cell(1,length(chs));
    max_digital = cellfun(@(x) strcmp(x,'max_digital'),aux_info(:,1));
    max_analog = cellfun(@(x) strcmp(x,'max_analog'),aux_info(:,1));
    analog_unit = cellfun(@(x) strcmp(x,'analog_unit'),aux_info(:,1));
    smpfilter_n = cellfun(@(x) strcmp(x,'smpfilter'),aux_info(:,1));
    uniques_sr = unique([x{:,2}]);
    sr_data = cell(1, length(uniques_sr)); % selected using sr index
    
    for ci = 1:length(chs)
        handles.ch_data{ci}.ch = chs(ci);
        handles.ch_data{ci}.conversion =  aux_info{max_analog,ci+1}/aux_info{max_digital,ci+1};
        handles.ch_data{ci}.unit =  aux_info{analog_unit,ci+1};
        handles.ch_data{ci}.label = labels{ci,1};
        handles.ch_data{ci}.sri = find(uniques_sr==x{ci,2});
        handles.ch_data{ci}.parsr = ['f' num2str(x{ci,2})];
        handles.ch_data{ci}.smpfilter_n = aux_info{smpfilter_n,ci+1};
        handles.ch_data{ci}.smpfilter = smpfilter2str(aux_info{smpfilter_n,ci+1});
    end

    handles.labels=labels(:,1);
    
    freq_line = handles.par.freq_line;
    for si = 1:length(uniques_sr)
        sr = uniques_sr(si);
        sr_data{si}.sr = sr;
        n_s = 2^floor(log2(handles.par.ffftlength*sr));
        nsamples = n_s*handles.par.n_blocks;
        sr_data{si}.n_s = n_s;
        sr_data{si}.win = barthannwin(n_s);
        [~ ,fs] = pwelch(zeros(1,nsamples),sr_data{si}.win,0,[],sr_data{si}.sr,'onesided');
        sr_data{si}.fs = fs;
        sr_data{si}.t = linspace(0,nsamples/sr,nsamples);

        %calc notchs for this freq
        notchs = 1:min(handles.par.num_notchs,floor((sr/2*0.9)/freq_line));
        K = 1; Z = []; P = [];

        for i= notchs
            w = freq_line*i/(sr/2);
            bw = handles.par.notch_width/(sr/2);
            [b_notch,a_notch] = iirnotch(w,bw);
            [zi,pi,ki] = tf2zpk(b_notch,a_notch);
            K = K * ki;
            Z(end+1:end+2) = zi;
            P(end+1:end+2) = pi;
        end

        [S,G] = zp2sos(Z,P,K);       % Convert to SOS
        parsr = ['f' num2str(sr)];
        sr_data{si}.notch.S = S;
        sr_data{si}.notch.G = G;
        if handles.par.custom_filter.(parsr).enable
            wpass = [handles.par.custom_filter.(parsr).bp1*2/sr handles.par.custom_filter.(parsr).bp2*2/sr];
            fstop_h = handles.par.fstop_h;
            [orden_pass, Wnorm_pass] = ellipord(wpass,[handles.par.fstop_l*2/sr fstop_h*2/sr],handles.par.Rp,handles.par.Rs);
            [z_pass,p_pass,k_pass] = ellip(orden_pass,par.Rp,par.Rs,Wnorm_pass);
            k_pass = k_pass * K;
            z_pass = [z_pass Z];
            p_pass = [p_pass P];
            [s_pass,g_pass] = zp2sos(z_pass,p_pass,k_pass);
            sr_data{si}.custom_filter.S = s_pass;
            sr_data{si}.custom_filter.G = g_pass;
        end
        for c = find(cellfun(@(x) x.sri==si,handles.ch_data))
            handles.ch_data{c}.sr = sr;
            handles.ch_data{c}.nmax = nsamples;
            handles.ch_data{c}.nupdated = 0;
            handles.ch_data{c}.buffer = zeros(1,nsamples,'int16');
            handles.ch_data{c}.cont = zeros(1,nsamples,'double');
            handles.ch_data{c}.cont_filtered = zeros(1,nsamples,'double');
            handles.ch_data{c}.psd = ones(length(fs),1,'double');
            handles.ch_data{c}.psd_filtered = ones(length(fs),1,'double');
        end
    end
    set(handles.N_lb,'String','0');
    
    handles.sr_data = sr_data;
    
    %GUI params
    if ~ isappdata(handles.fix_ypower_z_cb,'f3000uV')
        ylim_el = {'ypower','ypower_z','yraw','yfiltered'};
        ylim_cv = {'fix_ypower_z_cb','fix_ypower_cb','fix_yraw_cb','fix_yfiltered_cb'};
        for parstr = fieldnames(handles.par.x_power_manual)'
            for unit = {'uV','mV'}
                field = [parstr{1} unit{1}];
                for cb = ylim_cv
                    setappdata(handles.(cb{1}),field,false);
                end
                for et = ylim_el
                    setappdata(handles.([et{1} '_min']),field,0);
                    setappdata(handles.([et{1} '_max']),field,1);
                end
            end
        end
    end
    
    get_period = round(handles.par.ffftlength*0.9,3); %time with less precision
    plot_period = round((handles.par.n_blocks*get_period  )*1.1,3);
    handles.timer_buffer_loop = timer('Name','buffer_timer','TimerFcn',{@buffer_loop,hObject},'Period',plot_period,'ExecutionMode','fixedDelay');
    handles.timer_calc_loop = timer('Name','calc_timer','TimerFcn',{@calc_loop,hObject},'Period',get_period,'ExecutionMode','fixedSpacing');

    guidata(hObject, handles);
    start(handles.timer_buffer_loop)
    start(handles.timer_calc_loop)
    display_channel(1,hObject)
end



% --- Executes on button press in prev_ch.
function change_ch_Callback(hObject, inc)
% hObject    handle to prev_ch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    handles = guidata(hObject);
    ch_num = handles.channels_lb.Value;
    new_n = ch_num + inc;
    maxch = size(handles.channels_lb.String,1);
    if new_n > maxch
        new_n = 1;
    elseif new_n == 0
        new_n = maxch;
    end
    display_channel(new_n,hObject)
end

function save_b_Callback(hObject)
    h = guidata(hObject);
    ch_num = h.channels_lb.Value;
    selpath = uigetdir('choose folder');
    if h.save_all_cb.Value ==0
        saveas(h.cbmex_lfp,fullfile(selpath,[h.ch_data{ch_num}.label '.png']));
    else
        maxch = size(h.channels_lb.String,1);
        for ch = circshift(1:maxch,-4)
            display_channel(ch,hObject)
            drawnow
            saveas(h.cbmex_lfp,fullfile(selpath,[h.ch_data{ch}.label '.png']));
        end
    end
end
        
function channels_lb_Callback(hObject,eventdata,handles)
    display_channel(eventdata.Source.Value,hObject)
end

function string = smpfilter2str(n)
switch n
    case 0
        string = 'None';
    case 1
        string = '750Hz High pass';
    case 2
        string = '250Hz High pass';
    case 3
        string = '100Hz High pass';
    case 4
        string = '50Hz Low pass';
    case 5
        string = '125Hz Low pass';
    case 6
        string = '250Hz Low pass';
    case 7
        string = '500Hz Low pass';
    case 8
        string = '150Hz Low pass';
    case 9
        string = '10Hz-250Hz Band pass';
    case 10
        string = '2.5kHz Low pass';
    case 11
        string = '2kHz Low pass';
    case 12
        string = '250Hz-5kHz Band pass';
    otherwise
        string='unknown';
end
end


% --- Executes during object creation, after setting all properties.
function ax_CreateFcn(hObject)
    a=axtoolbar(hObject,{'export','datacursor','pan'	,'zoomin','zoomout','restoreview'});
    a.Visible = 'on';
end


% --- Executes on button press in fix_ypower_cb.
function fixscale_cb_Callback(hObject, eventdata)
% hObject    handle to fix_ypower_cb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    h = guidata(hObject);
    ch = h.channels_lb.Value;
    field = [h.ch_data{ch}.parsr  h.ch_data{ch}.unit];
    setappdata(hObject,field,hObject.Value);
    enable_editlines(hObject.Tag,hObject.Value,h)
    switch hObject.Tag
        case 'fix_ypower_cb'
            ax = h.spectrum;
            ymin = getappdata(h.ypower_min,field);
            ymax = getappdata(h.ypower_max,field);
        case 'fix_ypower_z_cb'
            ax = h.spectrum_zoom;
            ymin = getappdata(h.ypower_z_min,field);
            ymax = getappdata(h.ypower_z_max,field);
        case 'fix_yraw_cb'
            ax = h.time_raw;
            ymin = getappdata(h.yraw_min,field);
            ymax = getappdata(h.yraw_max,field);
        case 'fix_yfiltered_cb'
            ax = h.time_filtered;
            ymin = getappdata(h.yfiltered_min,field);
            ymax = getappdata(h.yfiltered_max,field);
    end
    if hObject.Value==1
        ylim(ax,[ymin ymax])
    else
        ylim(ax,'auto')
    end
end

function enable_editlines(checkbox_name,value,h)
    if value==1
        enable = 'on';
    else
        enable = 'off';
    end
    switch checkbox_name
        case 'fix_ypower_cb'
            set(h.ypower_min,'Enable',enable);
            set(h.ypower_max,'Enable',enable);
        case 'fix_ypower_z_cb'
            set(h.ypower_z_min,'Enable',enable);
            set(h.ypower_z_max,'Enable',enable);
        case 'fix_yraw_cb'
            set(h.yraw_min,'Enable',enable);
            set(h.yraw_max,'Enable',enable);
        case 'fix_yfiltered_cb'
            set(h.yfiltered_min,'Enable',enable);
            set(h.yfiltered_max,'Enable',enable);
    end
end

function setlim_Callback(hObject, e, h)
    new_value = str2num(e.Source.String);
    ch_num = h.channels_lb.Value;
	pars = h.ch_data{ch_num}.parsr;
	edit_tag = hObject.Tag;
    gui_par = [h.ch_data{ch_num}.parsr  h.ch_data{ch_num}.unit];

    if strcmp(edit_tag,'xpower_min') || strcmp(edit_tag,'xpower_max')
        if strcmp(edit_tag,'xpower_min')
            h.par.x_power_manual.(pars).min = new_value;
        else
            h.par.x_power_manual.(pars).max = new_value;
        end
        xlim(h.spectrum,[h.par.x_power_manual.(pars).min,h.par.x_power_manual.(pars).max]);
    end
    
    if strcmp(edit_tag,'xpower_z_min') || strcmp(edit_tag,'xpower_z_max') 
        if strcmp(edit_tag,'xpower_z_min')
            h.par.x_power_manual.(pars).min_zoom = new_value;
        else
            h.par.x_power_manual.(pars).max_zoom = new_value;
        end
        xlim(h.spectrum_zoom,[h.par.x_power_manual.(pars).min_zoom,h.par.x_power_manual.(pars).max_zoom]);
    end
    
    if strcmp(edit_tag,'ypower_min') || strcmp(edit_tag,'ypower_max') 
        if strcmp(edit_tag,'ypower_min')
            setappdata(h.ypower_min,gui_par,new_value);
        else
            setappdata(h.ypower_max,gui_par,new_value);
        end
        ylim(h.spectrum,[getappdata(h.ypower_min,gui_par),getappdata(h.ypower_max,gui_par)]);
    end
    
    if strcmp(edit_tag,'ypower_z_min') || strcmp(edit_tag,'ypower_z_max') 
        if strcmp(edit_tag,'ypower_z_min')
            setappdata(h.ypower_z_min,gui_par,new_value);
        else
            setappdata(h.ypower_z_max,gui_par,new_value);
        end
        ylim(h.spectrum_zoom,[getappdata(h.ypower_z_min,gui_par),getappdata(h.ypower_z_max,gui_par)]);
    end
    
    if strcmp(edit_tag,'yraw_min') || strcmp(edit_tag,'yraw_max') 
        if strcmp(edit_tag,'yraw_min')
            setappdata(h.yraw_min,gui_par,new_value);
        else
            setappdata(h.yraw_max,gui_par,new_value);
        end
        ylim(h.time_raw,[getappdata(h.yraw_min,gui_par),getappdata(h.yraw_max,gui_par)]);
    end
    
    if strcmp(edit_tag,'yfiltered_min') || strcmp(edit_tag,'yfiltered_max') 
        if strcmp(edit_tag,'yfiltered_min')
            setappdata(h.yfiltered_min,gui_par,new_value);
        else
            setappdata(h.yfiltered_max,gui_par,new_value);
        end
        ylim(h.time_filtered,[getappdata(h.yfiltered_min,gui_par),getappdata(h.yfiltered_max,gui_par)]);
    end
    
    
    
    guidata(hObject, h);
end
