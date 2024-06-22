function [] = makeGridPlot(beta0, N, tUB, sUB, tLB, sLB)
%MAKEGRIDPLOT Summary of this function goes here
close all;
load clown;

% all lambda values have cont steady state.
p = makeP(beta0, 1e-2, 'gridPlot');

if (tLB > 0) 
    p.usrlam = [tLB];
    p.nc.ilam=1; p.sol.ds=0.03; p.nc.dsmax=0.06; p.nc.lammin=0.0; p.sw.bifcheck=2; 
    p.nc.dlammax = 0.03; p.nc.dsmin=1e-10;
    p.nc.lammax=tLB * 1.1;
    p = cont(p);

    usrlam_points = p.branch(1,ismember(p.branch(4,:),p.usrlam));
    if isempty(usrlam_points)
        fprintf("Could not find starting point");
        return;
    end
    p = loadp('gridPlot',sprintf('pt%d',usrlam_points(1)));
end

if (sLB > 0) 
    p.usrlam = [sLB];
    p.nc.ilam=2; p.sol.ds=0.03; p.nc.dsmax=0.06; p.nc.lammin=0.0; p.sw.bifcheck=2; 
    p.nc.dlammax = 0.03; p.nc.dsmin=1e-10;
    p.nc.lammax=sLB * 1.1;
    p = cont(p);

    usrlam_points = p.branch(1,ismember(p.branch(4,:),p.usrlam));
    if isempty(usrlam_points)
        fprintf("Could not find starting point");
        return;
    end
    p = loadp('gridPlot',sprintf('pt%d',usrlam_points(1)));
end

thetaRange = linspace(tLB, tUB, N);
dTheta = thetaRange(2) - thetaRange(1);
sRange = linspace(sLB,sUB,N);
dS = sRange(2) - sRange(1);

x = getpte(p).';
dx = x(2)-x(1);

pu = cell(N);
puttn = zeros(N);
pustn = zeros(N);
agreement = zeros(N);

pu{1,1} = p.u;
%puttn(1, 1) = u_thetaTwoNorm(p ,1);
%pustn(1, 1) = u_thetaTwoNorm(p, 2);

dThetaMin = dTheta / 10;
dSMin = dS / 10;

for sIndex=1:N
    for tIndex=1:N
        if tIndex > 1 && ~isempty(pu{tIndex-1, sIndex})
            p.u = pu{tIndex-1, sIndex};
            p.u(p.nu+2) = sRange(sIndex);
            p.u(p.nu+1) = thetaRange(tIndex-1);
            p.u = nloop(p,p.u);

            p.usrlam = [thetaRange(tIndex)];
            p.nc.ilam=1; p.sol.ds=dThetaMin/10; p.nc.dsmax=dThetaMin; p.nc.lammin=thetaRange(tIndex)-dTheta; p.sw.bifcheck=2; 
            p.nc.dlammax = dThetaMin; p.nc.dsmin=1e-10;
            p.nc.lammax=thetaRange(tIndex)+dThetaMin*2;
            p.myParams.noConv = false;
            p = cont(p);

            usrlam_points = p.branch(1,ismember(p.branch(4,:),p.usrlam));
            if isempty(usrlam_points) || p.myParams.noConv
                fprintf("Fold at %d %d\n", tIndex, sIndex);
                continue;
            end
            p = loadp('gridPlot',sprintf('pt%d',usrlam_points(1)));
            pu{tIndex, sIndex} = p.u;

            %puttn(tIndex, sIndex) = u_thetaTwoNorm(p, 1);
            %pustn(tIndex, sIndex) = u_thetaTwoNorm(p, 2);
            fprintf("found %d %d from the left\n", tIndex, sIndex);
        end
        
        if sIndex > 1
            fprintf("n");
        end
        
        if sIndex > 1 && ~isempty(pu{tIndex, sIndex-1})
            p.u = pu{tIndex, sIndex-1};
            p.u(p.nu+2) = sRange(sIndex-1);
            p.u(p.nu+1) = thetaRange(tIndex);
            p.u = nloop(p,p.u);

            p.usrlam = [sRange(sIndex)];
            p.nc.ilam=2; p.sol.ds=dSMin/10; p.nc.dsmax=dSMin; p.nc.lammin=sRange(sIndex-1)-dSMin; p.sw.bifcheck=2; 
            p.nc.dlammax = dSMin; p.nc.dsmin=1e-10;
            p.nc.lammax=sRange(sIndex)+dSMin*2;
            p.myParams.noConv = false;
            p = cont(p,1000);

            usrlam_points = p.branch(1,ismember(p.branch(4,:),p.usrlam));
            if isempty(usrlam_points) || p.myParams.noConv
                fprintf("Fold at %d %d", tIndex, sIndex);
                continue;
            end
            p = loadp('gridPlot',sprintf('pt%d',usrlam_points(1)));
            
            if (~isempty(pu{tIndex, sIndex}))
                prev = pu{tIndex, sIndex};
                agreement(tIndex,sIndex) = norm(dx*(prev(1:p.nu) - p.u(1:p.nu))./ p.u(1:p.nu));
                if agreement(tIndex,sIndex) > 1e-3
                    figure(23);
                    plot(prev(1:p.nu) - p.u(1:p.nu));
                else
                    figure(24);
                    plot(prev(1:p.nu) - p.u(1:p.nu));
                end
            else
                pu{tIndex, sIndex} = p.u;
            end

            %puttn(tIndex, sIndex) = u_thetaTwoNorm(p, 1);
            %pustn(tIndex, sIndex) = u_thetaTwoNorm(p, 2);
            fprintf("found %d %d from the below\n", tIndex, sIndex);
            
        end
        
        figure(20);
        imagesc(thetaRange.', sRange.', puttn.')
        set(gca,'YDir','normal');
        colorbar;
        title("Norm of derivative of u with respect to theta 1");
        xlabel("theta 1");
        ylabel("theta 2");
        
        figure(22);
        imagesc(thetaRange, sRange,pustn.');
        set(gca,'YDir','normal');
        colorbar;
        title("Norm of derivative of u with respect to theta 2");
        xlabel("theta 1");
        ylabel("theta 2");
        
        
        figure(21);
        imagesc(thetaRange.', sRange.', agreement.')
        set(gca,'YDir','normal');
        title("Relative error between solutions ");
        colorbar;
        xlabel("theta 1");
        ylabel("theta 2");
        
        fprintf("Done with %d %d\n", tIndex, sIndex);
        
    end
end

end


function p = makeP(beta0, eps, fname)
    opts.fname = fname;
    
    opts.par.theta = 0.0;
    opts.par.alpha = 1;
    opts.par.d = 1/40;
    opts.par.epsilon = eps;
    opts.par.betaf = @(x,theta) beta0 + (1 - theta *cos(pi .* (x+0.5)))/25;
    opts.par.eta = @(x,theta) 1 - opts.par.betaf(x,theta);
    opts.L = 0.5;
    opts.Nx = 1e4;
    
    p = init_schk(opts);
    p.myParams.lScaling = 1;
    
    % Set continuation parameters
    
    p.nc.ilam=1;
    p.nc.dsmin = 1e-5;
    p.sol.ds=1e-3;
    p.nc.dsmax=1e-1;
    p.nc.dlammin = 1e-5;
    p.nc.dlammax = 1e-1;
    p.nc.lammin=0;
    p.nc.lammax=1;
    
    p.usrlam = [1];
    
    p.sw.eigssol = 0;
    p.sw.bifcheck = 1;
    p.sw.spcalc = 1;
    p.nc.tol = 1e-8;
    p.fsol.tol = 1e-10;
    p.fsol.fsol = 0;
    
    p.fsol.disp = 0;

end

