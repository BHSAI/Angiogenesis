function status = interuptFcn(t,y,flag,Z)
interupt_time = 60;
persistent INIT_TIME;
status = 0;
switch(flag)
    case 'init'
        INIT_TIME = tic;
    case 'done'
        clear INIT_TIME;
    otherwise
        elapsed_time = toc(INIT_TIME);
        if elapsed_time > interupt_time
            clear INIT_TIME;
            error('Stopped. Taking too long.');
            status=1;
        end
end
end