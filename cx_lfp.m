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
    N = getappdata(h.cbmex_lfp,'N');
    set(h.N_lb,'String',num2str(N));
    set(h.channels_lb,'Value',n);
    ch_info = getappdata(h.cbmex_lfp,'ch_info');
    par = getappdata(h.cbmex_lfp,'par');
    ch_data = getappdata(h.cbmex_lfp,'ch_data');
    if strcmp(ch_info(n).unit,'uV') && par.custom_filter.(ch_info(n).parsr).enable
        additional_curve = 'Notch and Custom Filter'; 
    else
        additional_curve = 'Notch'; 
    end
    log_psd = 10*log10(ch_data(n).psd);
    log_psd_filtered = 10*log10(ch_data(n).psd_filtered);
    set(h.channel_label,'String',sprintf('%s | ID: %d', ch_info(n).label,ch_info(n).ch))

    %set gui params
    set(h.info_channel,'String',{sprintf('Sampl. Rate: %1.0f Hz',ch_info(n).sr),...
        sprintf('Cont. Acq.: %s',ch_info(n).smpfilter)});
    
    pars = ch_info(n).parsr;
    set(h.xpower_min,'String',num2str(par.x_power_manual.(pars).min));
    set(h.xpower_max,'String',num2str(par.x_power_manual.(pars).max));
    set(h.xpower_z_min,'String',num2str(par.x_power_manual.(pars).min_zoom));
    set(h.xpower_z_max,'String',num2str(par.x_power_manual.(pars).max_zoom));  
    
    gui_par = [ch_info(n).parsr  ch_info(n).unit];
    
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
    sr_data = getappdata(h.cbmex_lfp,'sr_data');

    plot(h.spectrum,sr_data(ch_info(n).sri).fs,log_psd,'LineWidth',1.5);
    if getappdata(h.fix_ypower_cb,gui_par)
        yl = [getappdata(h.ypower_min,gui_par),getappdata(h.ypower_max,gui_par)];
    else
        ylim(h.spectrum,'auto');
        yl = ylim(h.spectrum);
        yl(1) = yl(1) - 0.05*(yl(2)-yl(1));  
    end
    plot(h.spectrum,sr_data(ch_info(n).sri).fs,log_psd_filtered,'r','LineWidth',0.75)
    ylim(h.spectrum,yl)
    h.spectrum.XMinorGrid = 'on'; h.spectrum.YMinorGrid = 'on';    
    xlim(h.spectrum,[par.x_power_manual.(pars).min,par.x_power_manual.(pars).max]);
    xlabel(h.spectrum,'Frequency (Hz)')
    ylabel(h.spectrum,'Power Spectrum (dB/Hz)')
    legend(h.spectrum,'Raw',additional_curve)
    
    %spectrum_zoom
    cla(h.spectrum_zoom); hold(h.spectrum_zoom,'on');
    plot(h.spectrum_zoom,sr_data(ch_info(n).sri).fs,log_psd,'LineWidth',1.5);
    if getappdata(h.fix_ypower_z_cb,gui_par)
        yl = [getappdata(h.ypower_z_min,gui_par),getappdata(h.ypower_z_max,gui_par)];
    else
        ylim(h.spectrum_zoom,'auto');
        yl = ylim(h.spectrum_zoom);
        yl(1) = yl(1) - 0.05*(yl(2)-yl(1));  
    end
    plot(h.spectrum_zoom,sr_data(ch_info(n).sri).fs,log_psd_filtered,'r','LineWidth',0.75)
    ylim(h.spectrum_zoom,yl)
    h.spectrum_zoom.XMinorGrid = 'on'; h.spectrum_zoom.YMinorGrid = 'on';
    xlim(h.spectrum_zoom,[par.x_power_manual.(pars).min_zoom,par.x_power_manual.(pars).max_zoom]);
    xlabel(h.spectrum_zoom,'Frequency (Hz)')
    ylabel(h.spectrum_zoom,'Power Spectrum (dB/Hz)')

    %time_raw
    cla(h.time_raw)
    plot(h.time_raw,sr_data(ch_info(n).sri).t,ch_data(n).cont,'LineWidth',1.2);
    h.time_raw.XMinorGrid = 'on'; h.time_raw.YMinorGrid = 'on';
    ylabel(h.time_raw,['Raw (' ch_info(n).unit ')'])
    xlabel(h.time_raw,'Time (sec)')
    xlim(h.time_raw,[sr_data(ch_info(n).sri).t(1) sr_data(ch_info(n).sri).t(end)])
    if getappdata(h.fix_yraw_cb,gui_par)
        ylim(h.time_raw,[getappdata(h.yraw_min,gui_par),getappdata(h.yraw_max,gui_par)])
    end
    %time_filtered
    cla(h.time_filtered)
    plot(h.time_filtered,sr_data(ch_info(n).sri).t,ch_data(n).cont_filtered,'r','LineWidth',1.2);
    h.time_raw.XMinorGrid = 'on'; h.time_filtered.YMinorGrid = 'on';
    ylabel(h.time_filtered,{additional_curve,['(' ch_info(n).unit ')']})
    xlabel(h.time_filtered,'Time (sec)')
    
    
    if getappdata(h.fix_yfiltered_cb,gui_par)
        ylim(h.time_filtered,[getappdata(h.yfiltered_min,gui_par),getappdata(h.yfiltered_max,gui_par)])
    end
    xlim(h.time_filtered,[sr_data(ch_info(n).sri).t(1) sr_data(ch_info(n).sri).t(end)])
    linkaxes([h.time_raw,h.time_filtered],'x')
end





function buffer_loop(this_timer,~,cbmex_lfp,handles)
%function that is called using a timer and update the gui data
    par = getappdata(cbmex_lfp,'par');
    stop(this_timer)
    [~, nsp_buffer] = cbmex('trialdata',1);
    if size(nsp_buffer,1) == 0
        pause(0.01)
        [~, nsp_buffer] = cbmex('trialdata',1);
    end
    keep_loading = false;
    buffer = getappdata(cbmex_lfp,'buffer');
    for ci = 1:length(par.channels)
        news = size(nsp_buffer{ci,3},1);
        if news == 102400 %by default max datapoints per channel 102400
            warning('lossing continuos data in nsp buffer')
        end
        to_update = min(news,buffer(ci).nmax*2-buffer(ci).nupdated);
        if to_update < news
            warning('lossing continuos data in matlab buffer')
        end
        buffer(ci).data(buffer(ci).nupdated+1:buffer(ci).nupdated+to_update) = nsp_buffer{ci,3}(1:to_update);
        buffer(ci).nupdated = to_update + buffer(ci).nupdated;
        
        
        if ci==1
            extra_data = mod(buffer(ci).nupdated,buffer(ci).nmax); %in the ideal case extra_Data should be 0
            next_call = par.ffftlength - extra_data/nsp_buffer{ci,2};
            next_call = round(next_call,3);
            this_timer.StartDelay = next_call;
            start(this_timer)
        end
        
        if buffer(ci).nmax > buffer(ci).nupdated
            keep_loading = true; %todos tienen que ser true!!!!!
        end
    end
    if keep_loading == false %all buffers are full
        sr_data = getappdata(cbmex_lfp,'sr_data');
        ch_info = getappdata(cbmex_lfp,'ch_info');
        ch_data = getappdata(cbmex_lfp,'ch_data');
        N = getappdata(cbmex_lfp,'N');
        N = N + 1; 
        bl_part = mod(N-1,par.n_blocks);
        for i = 1:length(par.channels)
            new_segment = (bl_part)*buffer(i).nmax+1:(bl_part+1)*buffer(i).nmax;
            ch_data(i).cont(new_segment) =  double(buffer(i).data(1:buffer(i).nmax))*ch_info(i).conversion;

            %copy extra samples to beggining of file and move nupdate
            buffer(i).data(1:buffer(i).nupdated - buffer(i).nmax) = buffer(i).data(buffer(i).nmax+1:buffer(i).nupdated);
            buffer(i).nupdated = buffer(i).nupdated - buffer(i).nmax;
            
            si = ch_info(i).sri ;  %sample rate
            psd = periodogram(ch_data(i).cont(new_segment(1:sr_data(si).fft_n_s)),sr_data(si).win,buffer(i).nmax,sr_data(si).sr,'onesided');
            ch_data(i).psd(:) = psd/N + ch_data(i).psd*(N-1)/N;

            if strcmp(ch_info(i).unit,'uV') && par.custom_filter.(ch_info(i).parsr).enable  
                s = sr_data(si).custom_filter.S;
                g = sr_data(si).custom_filter.G;
            else
                s = sr_data(si).notch.S;
                g = sr_data(si).notch.G;
            end
            ch_data(i).cont_filtered(new_segment)= filtfilt(s,g,ch_data(i).cont(new_segment));

            [psd ,~] = pwelch(ch_data(i).cont_filtered,sr_data(si).win,0,[],sr_data(si).sr,'onesided');
            ch_data(i).psd_filtered(:) = psd/N + ch_data(i).psd_filtered*(N-1)/N;       

        end
        setappdata(cbmex_lfp,'ch_data',ch_data);
        setappdata(cbmex_lfp,'N',N)
        if bl_part==(par.n_blocks-1) && handles.stop_refresh_cb.Value==false
            display_channel(handles.channels_lb.Value,cbmex_lfp)
        end
    end
    setappdata(cbmex_lfp,'buffer',buffer)
    
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
        stop(handles.timer_buffer);
        cbmex('close');
        delete(handles.timer_buffer);
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
    par = par_cb_lfp();
    if isempty(par.IP_address)
		[connection source] = cbmex('open');
	else
		[connection source] = cbmex('open','central-addr',par.IP_address, 'instance',0);
	end
    
    handles.new_data_loaded = false;
    chs = par.channels;
    pause(0.01)
    if ~isstr(chs)
        cbmex('mask',0,0)
        for j = chs
            pause(0.01)
            cbmex('mask',j,1)
        end
    else
        cbmex('mask',0,1)
    end
    pause(0.01)
    cbmex('trialconfig',1,'noevent');
    pause(0.01)
    [~, x] = cbmex('trialdata',1);
    while size(x,1) == 0
        pause(0.01)
        [~, x] = cbmex('trialdata',1);
    end
    chs = cell2mat(x(:,1));
    par.channels = chs;
    handles.cbmex_lfp.PaperPosition = handles.cbmex_lfp.Position;
    handles.cbmex_lfp.PaperPositionMode = 'auto';
    pause(0.01)
    labels = cbmex('chanlabel',chs);
    set(handles.channels_lb,'string',labels(:,1));
    
    aux_info = {};
    for c = chs(:)'
        pause(0.01)
        if isempty(aux_info)  
            aux_info = cbmex('config',c);
        else
            tmp = cbmex('config',c);
            aux_info = [aux_info, tmp(:,2)];
        end
    end
    
    max_digital = cellfun(@(x) strcmp(x,'max_digital'),aux_info(:,1));
    max_analog = cellfun(@(x) strcmp(x,'max_analog'),aux_info(:,1));
    analog_unit = cellfun(@(x) strcmp(x,'analog_unit'),aux_info(:,1));
    smpfilter_n = cellfun(@(x) strcmp(x,'smpfilter'),aux_info(:,1));
    uniques_sr = unique([x{:,2}]);
    ch_info = struct();
    for ci = 1:length(chs)
        ch_info(ci).ch = chs(ci);
        ch_info(ci).conversion =  aux_info{max_analog,ci+1}/aux_info{max_digital,ci+1};
        ch_info(ci).unit =  aux_info{analog_unit,ci+1};
        ch_info(ci).label = labels{ci,1};
        ch_info(ci).sri = find(uniques_sr==x{ci,2});
        ch_info(ci).sr = x{ci,2};
        ch_info(ci).parsr = ['f' num2str(x{ci,2})];
        ch_info(ci).smpfilter_n = aux_info{smpfilter_n,ci+1};
        ch_info(ci).smpfilter = smpfilter2str(aux_info{smpfilter_n,ci+1});
    end
    setappdata(handles.cbmex_lfp,'ch_info',ch_info);
    freq_line = par.freq_line;
    ch_data = struct();
    buffer = struct();
    sr_data = struct(); % selected using sr index
    
    
    for si = 1:length(uniques_sr)
        sr = uniques_sr(si);
        sr_data(si).sr = sr;
        n_s = ceil(par.ffftlength*sr); %longuitude for raw data
        fft_n_s = 2^floor(log2(par.ffftlength*sr)); 
        nsamples = n_s*par.n_blocks;
        sr_data(si).n_s = n_s;
        sr_data(si).fft_n_s = fft_n_s;
        sr_data(si).win = barthannwin(fft_n_s);
        [~ ,fs] = periodogram(zeros(1,fft_n_s),sr_data(si).win,fft_n_s,sr_data(si).sr,'onesided');
        sr_data(si).fs = fs;
        sr_data(si).t = linspace(0,nsamples/sr,nsamples);

        %calc notchs for this freq
        %notchs = 1:min(par.num_notchs,floor((sr/2*0.9)/freq_line));
        K = 1; Z = []; P = [];

        for i= 1:length(par.freq_line)
            w = par.freq_line(i)/(sr/2);
            %bw = par.notch_width/(sr/2);
            bw = w/par.notch_q;   
            [b_notch,a_notch] = iirnotch(w,bw);
            [zi,pi,ki] = tf2zpk(b_notch,a_notch);
            K = K * ki;
            Z(end+1:end+2) = zi;
            P(end+1:end+2) = pi;
        end

        [S,G] = zp2sos(Z,P,K);       % Convert to SOS
        parsr = ['f' num2str(sr)];
        sr_data(si).notch.S = S;
        sr_data(si).notch.G = G;
        
        if par.custom_filter.(parsr).enable
            wpass = [par.custom_filter.(parsr).bp1*2/sr par.custom_filter.(parsr).bp2*2/sr];
            fstop_h = par.fstop_h * par.custom_filter.(parsr).bp2;
            [orden_pass, Wnorm_pass] = ellipord(wpass,[par.fstop_l*2/sr fstop_h*2/sr],par.Rp,par.Rs);
            [z_pass,p_pass,k_pass] = ellip(orden_pass,par.Rp,par.Rs,Wnorm_pass);
            k_pass = k_pass * K;
            z_pass = [z_pass; Z'];
            p_pass = [p_pass; P'];
            [s_pass,g_pass] = zp2sos(z_pass,p_pass,k_pass);
            sr_data(si).custom_filter.S = s_pass;
            sr_data(si).custom_filter.G = g_pass;
        end
        for c = find(cellfun(@(x) x==si,{ch_info.sri}))
            buffer(c).nmax = n_s;
            buffer(c).nupdated = 0;
            buffer(c).data = zeros(1,2*n_s,'int16');
            ch_data(c).cont = zeros(1,nsamples,'double');
            ch_data(c).cont_filtered = zeros(1,nsamples,'double');
            ch_data(c).psd = ones(length(fs),1,'double');
            ch_data(c).psd_filtered = ones(length(fs),1,'double');
        end
    end
    set(handles.N_lb,'String','0');
    setappdata(handles.cbmex_lfp,'sr_data',sr_data);
    setappdata(handles.cbmex_lfp,'buffer',buffer);
    setappdata(handles.cbmex_lfp,'ch_data',ch_data);
    setappdata(handles.cbmex_lfp,'par',par);
    setappdata(handles.cbmex_lfp,'N',0);
    %GUI params
    if ~ isappdata(handles.fix_ypower_z_cb,'f3000uV')
        ylim_el = {'ypower','ypower_z','yraw','yfiltered'};
        ylim_cv = {'fix_ypower_z_cb','fix_ypower_cb','fix_yraw_cb','fix_yfiltered_cb'};
        for parstr = fieldnames(par.x_power_manual)'
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
    
    handles.timer_buffer = timer('Name','buffer_timer','Period',100,...
        'ExecutionMode','fixedSpacing','StartDelay',0.001);
    handles.timer_buffer.TimerFcn = {@buffer_loop,...
                                            handles.cbmex_lfp,handles};
    guidata(hObject, handles);
    display_channel(1,hObject)
    start(handles.timer_buffer)
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
    ch_info = getappdata(hObject,h.cbmex_lfp);
    ch_num = h.channels_lb.Value;
    selpath = uigetdir('choose folder');
    if h.save_all_cb.Value ==0
        saveas(h.cbmex_lfp,fullfile(selpath,[ch_info(ch_num).label '.png']));
    else
        maxch = size(h.channels_lb.String,1);
        for ch = circshift(1:maxch,-4)
            display_channel(ch,hObject)
            drawnow
            saveas(h.cbmex_lfp,fullfile(selpath,[ch_info(ch_num).label '.png']));
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
    ch_info = getappdata(h.cbmex_lfp,'ch_info');
    ch = h.channels_lb.Value;
    field = [ch_info(ch).parsr  ch_info(ch).unit];
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
    ch_info = getappdata(h.cbmex_lfp,'ch_info');
    par = getappdata(h.cbmex_lfp,'par');
	pars = ch_info(ch_num).parsr;
    
	edit_tag = hObject.Tag;
    gui_par = [pars  ch_info(ch_num).unit];

    if strcmp(edit_tag,'xpower_min') || strcmp(edit_tag,'xpower_max')
        if strcmp(edit_tag,'xpower_min')
            par.x_power_manual.(pars).min = new_value;
        else
            par.x_power_manual.(pars).max = new_value;
        end
        setappdata(hObject,h.cbmex_lfp,'par',par);
        xlim(h.spectrum,[par.x_power_manual.(pars).min,par.x_power_manual.(pars).max]);
    end
    
    if strcmp(edit_tag,'xpower_z_min') || strcmp(edit_tag,'xpower_z_max') 
        if strcmp(edit_tag,'xpower_z_min')
            par.x_power_manual.(pars).min_zoom = new_value;
        else
            par.x_power_manual.(pars).max_zoom = new_value;
        end
        setappdata(hObject,h.cbmex_lfp,'par',par);
        xlim(h.spectrum_zoom,[par.x_power_manual.(pars).min_zoom,par.x_power_manual.(pars).max_zoom]);
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
end
