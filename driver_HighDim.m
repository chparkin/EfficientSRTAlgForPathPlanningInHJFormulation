clear; %rng(1029);
% Solve the equation
%       u_t + 1_{x\neq xf}(x) (v(x)|grad(u)|-1) = 0
%       u(x,0) = i_{xf}(x)
% using Alex's splitting method
%
% Here H(x,p) = 1_{x\neq xf}(x) * (v(x)|p|-1),  g(x) = i_{xf}(x)
%      

% spatial point you are solving at
x_target = -0.5*ones(100,1);
xf = 0.5*ones(100,1);
xo = zeros(size(x_target));
xo(1) = 0.1;

% velocity: v(x) = 0.75 + 0.25*sin(x1)*sin(x2)*sin(x3)*...
v = @(x,t) 1-0.9*exp(-norm(x-xo,2)^2);
grad_v = @(x,t) 1.8*exp(-norm(x-xo,2)^2)*(x-xo);

% time you are solving at, and how many steps
t = 15; dt1 = 0.1; 
s = 0:dt1:t;
  
  
% set parameters        
sig = 1; tau = 0.2/sig; theta = 1; max_iter = 30000; tol = 5e-4;
gd_steps = 1; gd_rate = 0.2;

% solve equation
timer = tic;
% [u,x,p,howManyIter] = SolveHJB_LaxWithVel(x_target,xf,t,s,sig,tau,theta,max_iter,tol,gd_steps,gd_rate);
[u,x,p,howManyIter] = SolveHJB_LaxWithVelHighDim(x_target,xf,t,s,v,grad_v,sig,tau,theta,max_iter,tol,gd_steps,gd_rate);
TIME = toc(timer);
fprintf('Finished in %.2f seconds.\n',TIME);

%plot results 
% 
% xg = -5:0.05:5; yg = xg; 
% [X,Y] = ndgrid(xg,yg);
% figure(141);
% m=1;
% N = length(s)-1;
% circX = cos(0:2*pi/100:2*pi);
% circY = sin(0:2*pi/100:2*pi);
% for n = N+1:-5:1
% clf;hold on;
% contourf(X,Y,1-exp(-10*(X.^2 + (Y+0.2).^2)),50,'edgecolor','none'); 
% for m = 1:size(xC,2)
%     fill(xC(1,m)+r(m)*circX,xC(2,m)+r(m)*circY,[0.5,0,0],'edgecolor','none');
% end
% %C = cos(0:pi/200:2*pi); S = sin(0:pi/200:2*pi)-sin(1.8*s(n));
% %fill(C,S,'k');
% T = title(sprintf('$u = %.2f,\\,\\,t=%.2f$',u,t-s(n))); T.FontSize = 16; T.Interpreter = 'latex';
% plot(x(1,n:N+1),x(2,n:N+1),'k','linewidth',2);
% plot(xf(1),xf(2),'r.','markersize',20);
% plot(x_target(1),x_target(2),'g.','markersize',20);
% axis([-1.2 1.2 -1.2 1.2]);
% % axis off;
% pause(0.2);
% m=m+5;
% end
%%
xg = -5:0.05:5; yg = xg; 
[X,Y] = ndgrid(xg,yg);
figure(1123);clf; hold on;
contourf(X,Y,1-0.9*exp(-(X-0.1).^2 - Y.^2),50,'edgecolor','none');
plot(x(1,:),x(2,:),'k','linewidth',2);
plot(xf(1),xf(2),'r.','markersize',20);
plot(x_target(1),x_target(2),'.','markersize',20,'color',[0 0.7 0]); 
xlabel('$x_1$'); ylabel('$x_2$');
xticks(-1:0.5:1); yticks(-1:0.5:1);
ax = gca; ax.FontSize = 18; ax.TickLabelInterpreter = 'latex';
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex'; axis square;
% ax.YLabel.Rotation = 0;
axis([-1.5 1.5 -1.5 1.5]);
%%
figure(1124);clf;hold on; 
plot(s,fliplr(x(1,:)),'color','b','linewidth',2); 
plot(s,fliplr(x(2,:)),'color',rand(3,1),'linewidth',1.3); 
for i = 3:size(x,1)
   plot(s,fliplr(x(i,:)),'color',rand(3,1),'handlevisibility','off','linewidth',1.3); 
end
plot(s,fliplr(x(1,:)),'color','b','linewidth',2); 
plot(s,fliplr(x(2,:)),'color',rand(3,1),'linewidth',1.3);
plot(s,fliplr(x(1,:)),'color','b','linewidth',2,'handlevisibility','off');
L = legend({'$x_1(t)$','$x_i(t)$'}); L.FontSize = 18; L.Interpreter = 'latex'; L.Location = 'southeast';
xticks(0:5:15); yticks(-1:0.5:0.5);
xlabel('$t$'); ylabel('$x_i$');
ax = gca; ax.FontSize = 18; ax.TickLabelInterpreter = 'latex';
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.Rotation = 0;
axis([0 15 -1.4 0.7]);
