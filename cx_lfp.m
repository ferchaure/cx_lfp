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

% Last Modified by GUIDE v2.5 15-Jan-2016 21:25:40

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


% --- Executes just before cx_lfp is made visible.
function cx_lfp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cx_lfp (see VARARGIN)

% Choose default command line output for cx_lfp
handles.output = hObject;
% Update handles structure
setappdata(handles.axes1,'Ymanual',false);
guidata(hObject, handles);
start_adq(hObject)

% --- Outputs from this function are returned to the command line.
function varargout = cx_lfp_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);


function calc_loop(obj, event, hObject)
    handles = guidata(hObject);
    if handles.new_data_flag
        
        par = handles.par;
        handles.N = handles.N + 1; 
        N = handles.N;

        win = handles.win;
        for i = 1:length(par.channels)
            si = handles.sr4ch(i);  %sample rate
        	[psd ,~] = periodogram(handles.x{i},win{si},handles.data{si}.n_s,handles.data{si}.sr,'oneside');
            handles.PSD{i} = psd(1:length(handles.data{si}.fs))' /N + handles.PSD{i}*(N-1)/N;
        end
        handles.new_data_flag = false;
        guidata(hObject,handles);
        
        si = handles.sr4ch(handles.n_show);  %sample rate
        log_psd = log10(handles.PSD{handles.n_show})';
        plot(handles.axes1,handles.data{si}.fs,10*log_psd);
        Ymanual = getappdata(handles.axes1,'Ymanual');
        if Ymanual
           ylim(handles.axes1,[str2num(get(handles.ymin_et,'String')),str2num(get(handles.ymax_et,'String'))]) 
        end
        %power fit
        power_fit = log_psd(handles.data{si}.index_fit_1);
        p_all1 = polyfit(handles.data{si}.freqs_fit_1,power_fit,1);
        k1=p_all1(1);
        loga=p_all1(2);
        model = k1*handles.data{si}.freqs_fit_1+loga;
        [r2_all1 ~] = rsquare(power_fit,model);
         
        power_fit = log_psd(handles.data{si}.index_fit_2);
        p_all2 = polyfit(handles.data{si}.freqs_fit_2,power_fit,1);
        k2 = p_all2(1);
        loga = p_all2(2);
        model = k2*handles.data{si}.freqs_fit_2+loga;
        [r2_all2 ~] = rsquare(power_fit,model); 
        
        fit_text = sprintf('POWER FIT (CH %d) \n',par.channels(handles.n_show));
        fit_text = sprintf('%s %1.1f-%1.1f Hz: \n a = %1.2f (r^2 = %1.2f)\n',fit_text , par.fini_fit,par.fmid_fit,k1,r2_all1);
        fit_text = sprintf('%s %1.1f-%1.1f Hz: \n a = %1.2f (r^2 = %1.2f)\n',fit_text , par.fmid_fit,par.fmax_update,k2,r2_all2);
        set(handles.fit_text,'String',fit_text);
        
        set(handles.N_lb,'String',num2str(handles.N));
        xlim(handles.axes1,[par.fmin_disp,par.fmax_disp]);
        drawnow
        
    end



function buffer_loop(obj, event,hObject)
    handles = guidata(hObject);

    if handles.new_data_flag == false
        [t_ini_buffer x] = cbmex('trialdata',1);
        enogh_data = true;
        for i = 1:max(handles.sr4ch)
            ch = find(handles.sr4ch==i,1);
            if handles.data{i}.n_s > size(x{ch,3},1)
                enogh_data = false;
            end
        end
        
        if enogh_data
            for i = 1:size(x,1)
                si = handles.sr4ch(i);  %sample rate index
                handles.x{i} = x{i,3}(1:handles.data{si}.n_s);
            end
            handles.new_data_flag = true;
            guidata(hObject,handles);
        end
        
    else
        cbmex('trialdata',1);
    end

    


% --- Executes on button press in restart_button.
function restart_button_Callback(hObject, eventdata, handles)
% hObject    handle to restart_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    handles = guidata(hObject);
    handles.N = 0;
    guidata(hObject,handles);


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


function channel_et_Callback(hObject, eventdata, handles)
% hObject    handle to channel_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of channel_et as text
%        str2double(get(hObject,'String')) returns contents of channel_et as a double
handles = guidata(hObject);
ch_num = str2num(get(hObject,'String'));
chs = handles.par.channels;
if ~any(handles.par.channels == ch_num)
    set(hObject,'String',num2str(chs(handles.n_show)))
else
    handles.n_show = find(handles.par.channels == ch_num);
    set(handles.info_label,'String',handles.labels{handles.n_show});
    guidata(hObject,handles);
end

function fmax_et_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
handles.par.fmax_disp = str2num(get(hObject,'String'));
guidata(hObject,handles);

function fmin_et_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
handles.par.fmin_disp = str2num(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function et_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fmin_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close cx_lfp.
function cx_lfp_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to cx_lfp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
stop_adq(hObject)
delete(hObject);

function stop_adq(hObject)
    handles = guidata(hObject);
    stop(handles.timer_buffer_loop);
    stop(handles.timer_calc_loop);
    cbmex('close');
    delete(handles.timer_calc_loop);
    delete(handles.timer_buffer_loop);

    guidata(hObject, handles);

    
    
function start_adq(hObject)
    if exist('C:\Program Files\Blackrock Microsystems\NeuroPort Windows Suite', 'dir')
        addpath('C:\Program Files\Blackrock Microsystems\NeuroPort Windows Suite')
    else
        warning('Using cbmex simulator')
        addpath([pwd filesep 'sim'])
    end
    
    [connection source] = cbmex('open');
    handles = guidata(hObject);
    handles.par = par_cb_lfp();
    par = handles.par;
    handles.new_data_flag = false;
    handles.N = 0;
    handles.n_show = 1 ;
    chs = handles.par.channels;
    labels = cbmex('chanlabel',chs);
    handles.labels={labels{:,1}};
    
    set(handles.channel_et,'String',num2str(chs(handles.n_show)))
    set(handles.fmin_et,'String',num2str(handles.par.fmin_disp))
    set(handles.fmax_et,'String',num2str(handles.par.fmax_disp))
    set(handles.info_label,'String',handles.labels{handles.n_show});
        
    cbmex('mask',0,0)
    for j = chs
        cbmex('mask',j,1)
    end
    cbmex('trialconfig',1,'noevent');
    [t_ini_buffer x] = cbmex('trialdata',1);
    while size(x,1) == 0
        [t_ini_buffer x] = cbmex('trialdata',1);
    end
    handles.x = cell(1,length(chs));
    handles.PSD = cell(1,length(chs));
    uniques_sr = unique([x{:,2}]);
    data = cell(1, length(uniques_sr));
    handles.sr4ch = zeros(1,length(chs)); % using index of ch and sr 
    handles.win = cell(1,length(uniques_sr));
    for c = 1:length(uniques_sr)
        sr = uniques_sr(c);
        data{c}.sr = sr;
        
        n_s = 2^floor(log2(par.chunk_time*sr));
        data{c}.n_s = n_s;
        
        handles.sr4ch([x{:,2}]==sr) = c;
        handles.win{c} = barthannwin(n_s);
        [psd ,fs] = periodogram(zeros(1,n_s),handles.win{c},n_s,sr,'oneside');
        fs = fs(fs <= par.fmax_update);
        data{c}.fs = fs;
        
        data{c}.index_fit_1 = find(fs>par.fini_fit,1):find(fs>=par.fmid_fit,1);
        data{c}.freqs_fit_1 = log10(fs(data{c}.index_fit_1));

        data{c}.index_fit_2  = find(fs>par.fmid_fit,1):length(fs);
        data{c}.freqs_fit_2 = log10(fs(data{c}.index_fit_2));
        for i = find(handles.sr4ch==c)
            handles.x{i} = zeros(1,n_s);
            handles.PSD{i} = zeros(1,length(fs));
        end
    end

    handles.data = data;
    period = ceil(handles.par.chunk_time *100)/100; %time with less precision
    
    % UIWAIT makes cx_lfp wait for user response (see UIRESUME)
    % uiwait(handles.cx_lfp);
    
    handles.timer_buffer_loop = timer('Name','buffer_timer','TimerFcn',{@buffer_loop,hObject},'Period',period*1.1,'ExecutionMode','fixedDelay');
    handles.timer_calc_loop = timer('Name','calc_timer','TimerFcn',{@calc_loop,hObject},'Period',period*0.9,'ExecutionMode','fixedSpacing');

    guidata(hObject, handles);
    start(handles.timer_buffer_loop)
    start(handles.timer_calc_loop)

    
    
    
function [r2 rmse] = rsquare(y,f,varargin)
% Compute coefficient of determination of data fit model and RMSE
%
% [r2 rmse] = rsquare(y,f)
% [r2 rmse] = rsquare(y,f,c)
%
% RSQUARE computes the coefficient of determination (R-square) value from
% actual data Y and model data F. The code uses a general version of 
% R-square, based on comparing the variability of the estimation errors 
% with the variability of the original values. RSQUARE also outputs the
% root mean squared error (RMSE) for the user's convenience.
%
% Note: RSQUARE ignores comparisons involving NaN values.
% 
% INPUTS
%   Y       : Actual data
%   F       : Model fit
%
% OPTION
%   C       : Constant term in model
%             R-square may be a questionable measure of fit when no
%             constant term is included in the model.
%   [DEFAULT] TRUE : Use traditional R-square computation
%            FALSE : Uses alternate R-square computation for model
%                    without constant term [R2 = 1 - NORM(Y-F)/NORM(Y)]
%
% OUTPUT 
%   R2      : Coefficient of determination
%   RMSE    : Root mean squared error
%
% EXAMPLE
%   x = 0:0.1:10;
%   y = 2.*x + 1 + randn(size(x));
%   p = polyfit(x,y,1);
%   f = polyval(p,x);
%   [r2 rmse] = rsquare(y,f);
%   figure; plot(x,y,'b-');
%   hold on; plot(x,f,'r-');
%   title(strcat(['R2 = ' num2str(r2) '; RMSE = ' num2str(rmse)]))
%   
% Jered R Wells
% 11/17/11
% jered [dot] wells [at] duke [dot] edu
%
% v1.2 (02/14/2012)
%
% Thanks to John D'Errico for useful comments and insight which has helped
% to improve this code. His code POLYFITN was consulted in the inclusion of
% the C-option (REF. File ID: #34765).

if isempty(varargin); c = true; 
elseif length(varargin)>1; error 'Too many input arguments';
elseif ~islogical(varargin{1}); error 'C must be logical (TRUE||FALSE)'
else c = varargin{1}; 
end

% Compare inputs
if ~all(size(y)==size(f)); error 'Y and F must be the same size'; end

% Check for NaN
tmp = ~or(isnan(y),isnan(f));
y = y(tmp);
f = f(tmp);

if c; r2 = max(0,1 - sum((y(:)-f(:)).^2)/sum((y(:)-mean(y(:))).^2));
else r2 = 1 - sum((y(:)-f(:)).^2)/sum((y(:)).^2);
    if r2<0
    % http://web.maths.unsw.edu.au/~adelle/Garvan/Assays/GoodnessOfFit.html
        warning('Consider adding a constant term to your model') %#ok<WNTAG>
        r2 = 0;
    end
end

rmse = sqrt(mean((y(:) - f(:)).^2));



function ymax_et_Callback(hObject, eventdata, handles)
% hObject    handle to ymax_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ymax_et as text
%        str2double(get(hObject,'String')) returns contents of ymax_et as a double



function ymin_et_Callback(hObject, eventdata, handles)
% hObject    handle to ymin_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ymin_et as text
%        str2double(get(hObject,'String')) returns contents of ymin_et as a double


% --- Executes on button press in ymanual_tb.
function ymanual_tb_Callback(hObject, eventdata, handles)
% hObject    handle to ymanual_tb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ymanual_tb
if get(hObject,'Value')
    set(hObject, 'BackgroundColor',[0.992 0.612 0.557])
    set(handles.ymax_et,'Enable','on')
    set(handles.ymin_et,'Enable','on')
    set(handles.axes1,'YLimMode','manual')
    setappdata(handles.axes1,'Ymanual',true);
else
    set(hObject, 'BackgroundColor',[0.831 0.816 0.784])
    set(handles.ymax_et,'Enable','off')
    set(handles.ymin_et,'Enable','off')
    setappdata(handles.axes1,'Ymanual',false);
end

% --- Executes on button press in prev_ch.
function change_ch_Callback(hObject, eventdata, handles,inc)
% hObject    handle to prev_ch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
ch_num = str2num(get(hObject,'String'));
chs = handles.par.channels;
new_n = handles.n_show + inc;
if new_n > length(chs)
    new_n = 1;
elseif new_n == 0
    new_n = length(chs);
end
handles.n_show = new_n;
set(handles.info_label,'String',handles.labels{handles.n_show});
set(handles.channel_et,'String',num2str(chs(handles.n_show)))
guidata(hObject,handles);


function save_current_b_Callback(hObject, eventdata, handles)
    save_channel_spectrum(handles.n_show,handles)
    
function save_all_b_Callback(hObject, eventdata, handles)
    stop(handles.timer_buffer_loop);
    stop(handles.timer_calc_loop);
    save_channel_spectrum(1:length(handles.par.channels),handles)
    start(handles.timer_buffer_loop)
    start(handles.timer_calc_loop)
    
function save_channel_spectrum(ci,handles)       
        if ~exist(date,'dir')
            mkdir(date)
        end
        par = handles.par;
        handles.N = handles.N + 1; 
        N = handles.N;
        
        for ch = ci 
            fig = figure('Visible','Off') ; 
            si = handles.sr4ch(ch);  %sample rate
            log_psd = log10(handles.PSD{ch})';
            plot(handles.data{si}.fs,10*log_psd);
            Ymanual = getappdata(handles.axes1,'Ymanual');
            if Ymanual
               ylim([str2num(get(handles.ymin_et,'String')),str2num(get(handles.ymax_et,'String'))]) 
            end
            %power fit
            power_fit = log_psd(handles.data{si}.index_fit_1);
            p_all1 = polyfit(handles.data{si}.freqs_fit_1,power_fit,1);
            k1=p_all1(1);
            loga=p_all1(2);
            model = k1*handles.data{si}.freqs_fit_1+loga;
            [r2_all1 ~] = rsquare(power_fit,model);

            power_fit = log_psd(handles.data{si}.index_fit_2);
            p_all2 = polyfit(handles.data{si}.freqs_fit_2,power_fit,1);
            k2 = p_all2(1);
            loga = p_all2(2);
            model = k2*handles.data{si}.freqs_fit_2+loga;
            [r2_all2 ~] = rsquare(power_fit,model); 

            fit_text = sprintf('Ch %d:',par.channels(ch));
            fit_text = sprintf('%s %1.1f-%1.1fHz: a=%1.2f (r^2=%1.2f) -|-',fit_text , par.fini_fit,par.fmid_fit,k1,r2_all1);
            fit_text = sprintf('%s %1.1f-%1.1fHz: a=%1.2f (r^2=%1.2f)',fit_text , par.fmid_fit,par.fmax_update,k2,r2_all2);
            title(fit_text);

            xlim(handles.axes1,[par.fmin_disp,par.fmax_disp]);
            print(fig,[date filesep handles.labels{ch}],'-dpng')
        end
