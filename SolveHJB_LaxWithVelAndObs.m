function [u,x,p,howManyIter] = SolveHJB_LaxWithVelAndObs(x_target,xf,t,s,v,grad_v,xC,r,sig,tau,theta,max_iter,tol,gd_steps,gd_rate)
%
% This function runs algorithm 1 from the manuscript
%   to actually resolve and optimal path (and also
%   a value of the solution to the Hamilton-Jacobi
%   equation)
%
%
if nargin < 4
    sig = 1;
    tau = 0.2;
    theta = 1;
    max_iter = 1000;
    tol = 1e-6;
    gd_steps = 1;
    gd_rate = 0.1;
end

% Define some ambient parameters and initialize randomly
% with half the points near the starting point and half near the ending
% point
dim = size(xf,1);
N = length(s)-1;
x = xf + 0.1*(2*rand(dim,N+1)-1);
x(:,round((N-1)/2+1):N+1) = x(:,round((N-1)/2+1):N+1) -xf + x_target;
x(:,N+1) = x_target;
z = x;
p = 0.1*(2*rand(dim,N+1)-1);
p(:,1) = 0;

%obstacle parameters
aO = 100;
O = zeros(1,size(x,2)); % obstacle function values
gO = zeros(dim,size(x,2)); % gradient of the obstacle function


%approx indicator function of final point
a = 100;
chi = @(x) 1 - exp(-a*norm(x-xf,2)^2);
chix = @(x) 2*a*(x-xf)*exp(-a*norm(x-xf,2)^2);

% begin the interation
u = 100;
for k = 1:max_iter
    % store previous iterates
    x1 = x; p1 = p; u1 = u;
    
    % minimization over p vectors and resolution of obstacles
    for j = 2:N+1
        % resolve obstacles
        if isempty(r)
            O(j) = 1;
            gO(:,j) = [0;0];
        else
            [D,m] = min(sqrt((x(1,j)-xC(1,:)).^2 + (x(2,j)-xC(2,:)).^2)-r);
            O(j) = 1/2 + (1/2)*tanh(aO*D);
            gO(:,j) = ((aO/2)*(sech(aO*D)^2)/(D+r(m)+1e-8))*(x(:,j)-xC(:,m));
        end
        b = p(:,j) + sig*(z(:,j)-z(:,j-1));
        p(:,j) = max(0,1-(sig*O(j)*(s(j)-s(j-1))*chi(x(:,j))*v(x(1,j),x(2,j),s(j)))/norm(b,2))*b;
    end
    
    % minimization over x vectors
    x(:,1) = xf;
    for j = 2:N
        w = x(:,j) - tau*(p(:,j) - p(:,j+1));
        for l = 1:gd_steps
            x(:,j) = x(:,j) - gd_rate*(-(s(j+1)-s(j))*(chi(x(:,j))*norm(p(:,j),2)*(O(j)*grad_v(x(1,j),x(2,j),s(j))+v(x(1,j),x(2,j),s(j))*gO(:,j)) + (norm(p(:,j),2)*O(j)*v(x(1,j),x(2,j),s(j)) - 1)*chix(x(:,j))) + (1/tau)*(x(:,j) - w));
        end
        
    end
    
    % update z
    for j = 1:N+1
        z(:,j) = x(:,j) + theta*(x(:,j) - x1(:,j));
    end
    
    % see how much x and p changed
    change = max(max(abs(x(:)-x1(:))),max(abs(p(:)-p1(:))));
    
    %calculate u
    u = 0;
    for j = 2:N+1
        u = u + p(:,j)'*(x(:,j)-x(:,j-1))-(s(j)-s(j-1))*chi(x(:,j))*(O(j)*v(x(1,j),x(2,j),s(j))*norm(p(:,j),2)-1);
    end
    
    % see how much u changed (this was used to determine convergence 
    % in an earlier version, but now is simply a check to see how things
    % are going)
    uchange = abs(u-u1);
    
    % break if convergence tolerance is reached
    if  change < tol
        break;
    end
    
    % print progress
    if mod(k,1000)==0
        fprintf('Iteration %i complete, change in (x,p) = %.4e, change in u = %.4e\n',k,change,uchange);    
    end
    
    % halve the gradient descent rate each 1000 iterations after 5000
    if mod(k,1000)==0 && k >= 5000
        gd_rate = 0.5*gd_rate;
    end
    
end

%final computation of u
u = norm(x(:,1)-xf,2);
for j = 2:N+1
    u = u + p(:,j)'*(x(:,j)-x(:,j-1))-(s(j)-s(j-1))*chi(x(:,j))*(v(x(1,j),x(2,j),s(j))*norm(p(:,j),2)-1);
end
howManyIter = k;


end

