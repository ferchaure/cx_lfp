clear all
close all
clc
addpath('C:\Program Files\Blackrock Microsystems\NeuroPort Windows Suite')
[connection source] = cbmex('open');

hand1=figure(45654645);
chans = [2 105];
nchans = numel(chans);
% create mask for channels?
labels = cbmex('chanlabel',chans);
cbmex('mask',0,0)
for i=chans
    cbmex('mask',i,1)
end
%%
cbmex('trialconfig',1,'noevent');

t_read0 = tic;
read_buffer = 1;
slack=0.05;
while(ishandle(hand1))
    t_read=toc(t_read0);
    if t_read >=read_buffer+slack
        [t_ini_buffer data] = cbmex('trialdata',1);
        if read_buffer*data{1,2} > size(data{1,3},1)
            disp('menos muestras')
        else
            plot((1:read_buffer*data{1,2})/data{1,2},data{1,3}(1:read_buffer*data{1,2}))
            
            drawnow
        end
        t_read0 = tic;
    end
    pause(0.01)
end


cbmex('close');
% while(ishandle(hand1))
%     t_read=toc(t_read0);
%     if t_read >=read_buffer
%         [t_ini_buffer data] = cbmex('trialdata',1);
%         num_chan = size(data,1);
%         if (ishandle(hand1))    
%             for ichan=1:nchans
%                 clf(hand1)
%                 sr = data{ichan,2};
%                 x = data{ichan,3}(1:read_buffer*sr);
%                 subplot(2,1,1)
%                 plot((1:read_buffer*sr)/sr,x)
%                 title(sprintf('Channel %d raw',data{ichan,1}))
%             end
%         end
%         t_read0 = tic;
%     end 
% end
% cbmex('close');
