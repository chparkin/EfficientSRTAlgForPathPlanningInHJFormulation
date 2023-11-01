function [u,x,p,howManyIter] = SolveHJB_LaxWithVelHighDim(x_target,xf,t,s,v,grad_v,sig,tau,theta,max_iter,tol,gd_steps,gd_rate)
%SOLVEHJB Summary of this function goes here
%   Detailed explanation goes here
if nargin < 4
    sig = 1;
    tau = 0.2;
    theta = 1;
    max_iter = 1000;
    tol = 1e-6;
    gd_steps = 1;
    gd_rate = 0.1;
end
% XXX = 1;
xg = -5:0.05:5; yg = xg;
[X,Y] = ndgrid(xg,yg);

dim = size(xf,1);
N = length(s)-1;


x = xf + 0.1*(2*rand(dim,N+1)-1);
x(:,round((N-1)/2+1):N+1) = x(:,round((N-1)/2+1):N+1) -xf + x_target;
% x = 2*rand(2,N+1)-1;
x(:,N+1) = x_target;
z = x;
p = 0.1*(2*rand(dim,N+1)-1);
p(:,1) = 0;

a = 100;
chi = @(x) 1 - exp(-a*norm(x-xf,2)^2);
chix = @(x) 2*a*(x-xf)*exp(-a*norm(x-xf,2)^2);
%chix = @(x) 2*a*(x-xf)*((norm(x-xf,2)>0.2)*max(exp(-a*norm(x-xf,2)^2),1e-6) + (norm(x-xf,2)<=0.2)*exp(-a*norm(x-xf,2)^2));
% chi = @(x) 1*(norm(x-xf,2)>1/a) + a*norm(x-xf,2)*(norm(x-xf,2)<=1/a);
% chix = @(x) zeros(size(x))*(norm(x-xf,2)>1/a) + (a*(x-xf)/(1e-8+norm(x-xf,2)))*(norm(x-xf,2)<=1/a);
% chi = @(x) 1 - exp(-a*norm(x-xf,2));
% chix = @(x) (a*(x-xf)/(1e-8+norm(x-xf,2)))*exp(-a*norm(x-xf,2));
u = 100;
for k = 1:max_iter
    x1 = x; p1 = p; u1 = u;
    
    % minimization over p vectors
    for j = 2:N+1
        b = p(:,j) + sig*(z(:,j)-z(:,j-1));
        p(:,j) = max(0,1-(sig*(s(j)-s(j-1))*chi(x(:,j))*v(x(:,j),s(j)))/norm(b,2))*b;
    end
    % minimization over x vectors
    
    x(:,1) = xf;
    for j = 2:N
        w = x(:,j) - tau*(p(:,j) - p(:,j+1));
        for l = 1:gd_steps
            x(:,j) = x(:,j) - gd_rate*(-(s(j)-s(j-1))*(chi(x(:,j))*norm(p(:,j),2)*grad_v(x(:,j),s(j)) + (norm(p(:,j),2)*v(x(:,j),s(j)) - 1)*chix(x(:,j))) + (1/tau)*(x(:,j) - w));
        end

    end
    
    % update z
    for j = 1:N+1
        z(:,j) = x(:,j) + theta*(x(:,j) - x1(:,j));
    end
    
    change = max(max(abs(x(:)-x1(:))),max(abs(p(:)-p1(:))));
    
    %calculate u
    u = 0;
    for j = 2:N+1
        u = u + p(:,j)'*(x(:,j)-x(:,j-1))-(s(j)-s(j-1))*chi(x(:,j))*(v(x(:,j),s(j))*norm(p(:,j),2)-1);
    end
    
    uchange = abs(u-u1);
    
    if  change < tol
        break;
    end
    
%     figure(112);
%     if (mod(k,100)==1 && k < 5000) || (mod(k,500)==1 && k>5000)
%         clf;hold on;
%         contourf(X,Y,1+(9/10)*(sin(X).*sin(Y)),50,'edgecolor','none');
%         plot(x(1,:),x(2,:),'k','linewidth',2);
%         plot(xf(1),xf(2),'.','markersize',25,'color',[1 0 0]);
%         T = text(xf(1)-0.7,xf(2)-0.7,'$B$'); T.Interpreter = 'latex'; T.FontSize = 20; T.Color = [0.8 0 0];
%         plot(x_target(1),x_target(2),'.','markersize',25,'color',[0.6 0 0.8]);
%         T = text(x_target(1)-0.7,x_target(2)-0.7,'$A$'); T.Interpreter = 'latex'; T.FontSize = 20; T.Color = [0.6 0 0.8];
%         %C = cos(0:pi/200:2*pi); S = sin(0:pi/200:2*pi)-sin(1.8*s(n));
%         %fill(C,S,'k');
%         T = title(sprintf('Iteration: $%i$',k-1)); T.FontSize = 16; T.Interpreter = 'latex';
%         axis([-5 5 -5 5]);
%         axis off;
%         %         print(sprintf('PDHGDemo%i',XXX),'-dpng'); XXX = XXX+1;
%         pause(0.1);
%     end

    if mod(k,1000)==0
        fprintf('Iteration %i complete, change in (x,p) = %.4e, change in u = %.4e\n',k,change,uchange);
        if k >= 2000
            gd_rate = 0.5*gd_rate;
        end
    end
        
end

u = norm(x(:,1)-xf,2);
for j = 2:N+1
    u = u + p(:,j)'*(x(:,j)-x(:,j-1))-(s(j)-s(j-1))*chi(x(:,j))*(v(x(:,j),s(j))*norm(p(:,j),2)-1);
end
howManyIter = k;


end

