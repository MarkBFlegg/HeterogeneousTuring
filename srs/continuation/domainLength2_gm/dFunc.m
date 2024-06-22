function [status] = dFunc(t, y, flag)
    if flag == "init" || flag == "done"
        status = 0;
        return
    end
    dtnorm = 1;
    if dtnorm < 1e-9
        status = 1;
    else
        status = 0;
    end
    if t > 100 && dtnorm > 1
        status = 1;
    end
    figure(3)
    clf;
    hold on;
    plot(y);
    title(sprintf("time = %f, dt = %e ", t,dtnorm));
    drawnow;
    return
    
end

