% Solve the equation
%       u_t + 1_{x\neq xf}(x) (v(x)|grad(u)|-1) = 0
%       u(x,0) = i_{xf}(x)
% using method from manuscript
%
% Here H(x,p) = 1_{x\neq xf}(x) * (v(x)|p|-1),  g(x) = i_{xf}(x)
%
clear;
%% load which example you want to see
load ex1.mat
rng(SEED);

% To make your oen examples, you can modify the variables that are loaded
% from the files ex1.mat or ex2.mat. Parameters are:
%
%   x0 - initial point
%   xf - final point
%   t - travel time
%   v(x,y) - velocity function
%   grad_v(x,y) - gradient of velocity function (needs to be hard-coded)
%   dt - time discretization parameter
%   xO - coordinates for centers of obstacles - [x(1,m);x(2,m)] is the 
%               center of the mth obstacle
%   rO - radii of the obstacles, rO(m) is the radius of the mth obstacle
%   obs_tol - observation tolerance (how close the vehicle needs to be to 
%                                    spot an obstacle)
%
%   All other things included in the data files are variables used in the
%   solver or helper variables for plotting, and probably shouldn't be
%   modified.
%
%   To generate new obstacles, can use the function generateDisjointCircles
%
%   Changing the RNG seed will not meaningfully change results except in
%   cases where the optimal path is not unique (in which case the resolved
%   path depends on the initialization). 

%% perform optimization to find optimal path
timer2 = tic;
while norm(xPath(:,end)-xf,2)> end_of_path_tol
    
    % solve equation with knowledge of current obstacles
    timer = tic;
    [u(l),x{l},p,howManyIter] = SolveHJB_LaxWithVelAndObs(x_target,xf,t,s,v,grad_v,xCD,rD,sig,tau,theta,max_iter,tol,gd_steps,gd_rate);
    TIME = toc(timer);
    fprintf('Resolved path %i in %.2f seconds.\n',l,TIME);
    
    % travel until you hit an obstacle
    for j = size(x{l},2)-1:-1:1
        xPath(:,m) = x{l}(:,j);
        if ~isempty(rO)
        [D,e] = min(sqrt((x{l}(1,j)-xC(1,:)).^2 + (x{l}(2,j)-xC(2,:)).^2)-r);
        t = t-dt;
        if D < obs_tol
           for g = 1:size(xCO,2)
              if isequal(xCO(:,g),xC(:,e)) 
                  break;
              end
           end
           dInds(end+1) = g;
           dTimes(end+1) = (m-1)*dt;
           xCD(:,end+1) = xC(:,e);
           rD(:,end+1) = r(e);
           xC(:,e)=[]; r(e) = [];
           x_target = x{l}(:,j);
           l = l+1;
           fprintf('Discovered Obstacle at time %.2f. Re-computing path with t = %.2f, x_target = [%.2f;%.2f]\n',dTimes(end),t,x_target(1),x_target(2));
           m=m+1;
           break;
        end
        end
        m=m+1;
        if norm(xPath(:,end)-xf,2)<= end_of_path_tol
            break;
        end
    end
end
TIME2 = toc(timer2);
fprintf('Resolved full path in %.2f seconds. Discovered %i obstacle(s).\n',TIME2,length(dTimes));
%% plot results
figNum = 314159;
AX = [-1.2 1.2 -1.2 1.2];
PRINT = 0;
PRINTSTRING = 'SRTPlannerSplitDecision2_';
RecompInds = plotPath(u,x,xPath,dt,xCO,rO,dInds,dTimes,obs_tol,x0,xf,v,PRINT,PRINTSTRING,AX,figNum,vFlag);