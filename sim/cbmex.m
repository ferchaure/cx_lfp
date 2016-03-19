function varargout = cbmex(command, varargin)
persistent active_ch

max_ch = 180;
%ch second freq:
ch_other_freq = 100:102;

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
        varargout{2} = cell(num_act_ch, 3);
        for i = 1 : num_act_ch
            varargout{2}{i,1} = int32(channels(i));
            if ismember(channels(i),ch_other_freq)
                varargout{2}{i,2} = single(1000);
                cant_n = floor(2^15*1000/30000);
                varargout{2}{i,3} =  single(sin((1:cant_n)*i/100)'+rand(cant_n,1)); %example
            else
                varargout{2}{i,2} = single(30000);
                varargout{2}{i,3} = single(sin((1:2^15)*i/100)'+rand(2^15,1)); %example
            end
            
            
            
        end


    otherwise
        varargout{1} = 0;
end
end
