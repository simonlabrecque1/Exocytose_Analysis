function newEvents = removeRecEvents2(events,th_dist,th_time)

% Removes recurent events
% events  : array of events positions, time and intensity: x y t I
% th_dist : threshold in xy for two events to be considered separate
% th_time : threshold in time for two events to be considered separate

% tic
events = sortrows(events,-4);
newEvents = events(1,:);
it_new =1;
%h = waitbar(0,'please wait');
num = size(events,1);
for it = 2:num
%    waitbar(it./num);
    %[it num]
    timeDiff = abs(newEvents(:,3) -events(it,3)) < th_time;
    if (max(timeDiff)==1)
        xDiff = abs(newEvents(timeDiff,1) - events(it,1)) < th_dist;
        if (max(xDiff)==1)
            temp = newEvents(timeDiff,2);
            yDiff = abs(temp(xDiff) -events(it,2)) < th_dist; 
            if(max(yDiff)==1)
                %we don't take that point
            else
                it_new=it_new+1;
                newEvents(it_new,:) = events(it,:);   
            end
        else
            it_new=it_new+1;
            newEvents(it_new,:) = events(it,:);
        end
    else
        it_new=it_new+1;
        newEvents(it_new,:) = events(it,:);
    end
            
end
% toc
%close(h);
        
    