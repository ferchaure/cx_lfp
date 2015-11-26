function varargout = cbmex(command, varargin)
persistent active_ch
max_ch = 180;


if isempty(active_ch)
    active_ch = ones(max_ch,1);
end

switch command
    case 'open'
        varargout{1} = 0;
        varargout{2} = 0;
    case 'chanlabel'
        chs = varargin{1};
        varargout{1} = cell(length(chs),2);
        for i = 1:length(chs)
            varargout{1}{i,1} = ['sim_ch_' num2str(chs(i))];
        end
    case 'mask'
        if varargin{1}==0
            active_ch = zeros(max_ch,1) + varargin{2};
        else
            active_ch(varargin{1}) =  varargin{2};
        end    
    case 'trialdata'
        varargout{1} = 0;
        num_act_ch = nnz(active_ch);
        channels = find(active_ch);
        varargout{2} = cell(num_act_ch,3);
        for i =1 : num_act_ch
            varargout{2}{i,1} = channels(i);
            varargout{2}{i,2} = 30000;
            varargout{2}{i,3} = sin((1:2^15)*i/100)'+rand(2^15,1); %no estoy pensando como seria el espectro de esto... pero que son distintos seguro
        end
    otherwise
        varargout{1} = 0;
end
end
