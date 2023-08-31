function [RecompInds] = plotPath(u,x,xPath,dt,xCO,rO,dInds,dTimes,obs_tol,x0,xf,v,PRINT,PRINTNAME,AX,figNum,vFlag)

% a helper function to do all the plotting

if nargin <= 14
    figNum = 1;
    dx = abs(x0(1)-xf(1));
    dy = abs(x0(2)-xf(2));
    minx = min(x0(1),xf(1)) - dx/5;
    maxx = max(x0(1),xf(1)) + dx/5;
    miny = min(x0(2),xf(2)) - dy/5;
    maxy = max(x0(2),xf(2)) + dy/5;
    AX = [minx maxx miny maxy];
    vFlag = 0;
end
if nargin <= 12
    PRINT = 0;
    PRINTNAME = [];
end


RecompInds = [];

figure(figNum);
circX = cos(0:2*pi/100:2*pi);
circY = sin(0:2*pi/100:2*pi);
endTime = max(dTimes)+u(end);
l=1;
for n = 1:1:length(xPath)
    clf;hold on;
    axis(AX);
    if vFlag
        [X,Y] = ndgrid(AX(1):(AX(2)-AX(1))/100:AX(2),AX(3):(AX(4)-AX(3))/100:AX(4));
        contourf(X,Y,v(X,Y),50,'edgecolor','none');
    end
    for m = 1:size(xCO,2)
        COLOR = [0.5 0 0];
        DISC = find(dInds == m);
        if ~isempty(DISC)
            if (n-1)*dt>=dTimes(DISC)
                COLOR = [0 0 0.7];
            end
        end
        fill(xCO(1,m)+rO(m)*circX,xCO(2,m)+rO(m)*circY,COLOR,'edgecolor','none');
%          T = text(xCO(1,m),xCO(2,m),sprintf('%i',m)); T.FontSize = 15; T.Color = [0.9 0.9 0.9];
    end
    %C = cos(0:pi/200:2*pi); S = sin(0:pi/200:2*pi)-sin(1.8*s(n));
    %fill(C,S,'k');
    plot(xPath(1,n)+obs_tol*circX,xPath(2,n)+obs_tol*circY,'color',[0.5 0.5 0.5],'linewidth',1.5);
    P = plot(x{l}(1,:),x{l}(2,:),'--','linewidth',2,'color','m');
    plot(xPath(1,1:n),xPath(2,1:n),'k','linewidth',2);
    plot(xPath(1,n),xPath(2,n),'k.','markersize',20);
    plot(xf(1),xf(2),'r.','markersize',20);
    plot(x0(1),x0(2),'g.','markersize',20);
    for m = 1:length(dTimes)
        if (n-1)*dt>=dTimes(m)
            ind = round(dTimes(m)/dt)+1;
            plot(xPath(1,ind),xPath(2,ind),'c.','markersize',20);
        end
    end
    axis off;
     T = title(sprintf('$t=%.2f$\n',(n-1)*dt)); T.FontSize = 16; T.Interpreter = 'latex';
    if any(abs(dTimes-(n-1)*dt)<dt/3)
         T = title(sprintf('$t=%.2f$\n Discovered Obstacle - Recomputing Path',(n-1)*dt)); T.FontSize = 16; T.Interpreter = 'latex';
        RecompInds(end+1) = n-1;
        
        %%%% FOR PICS FROM PAPER
%         delete(P);
%         P = plot(x{l+1}(1,:),x{l+1}(2,:),'--','linewidth',2,'color','m');
%         plot(xPath(1,ind),xPath(2,ind),'c.','markersize',20);
%         print(sprintf('ex2_pic%i',l),'-dpng');
%         delete(P);
        %%%%%%%%%%%%%%%%%%%%%%
        
        pause(1);
        l=l+1;
        T = title(sprintf('$t=%.2f$\n New Path Computed',(n-1)*dt)); T.FontSize = 16; T.Interpreter = 'latex';
        pause(0.4);
        %         set(P,'Visible','off');
%         plot(x{l}(1,:),x{l}(2,:),'--','linewidth',2,'color','m');
%         pause(1);
    end
    if (n-1)*dt>=endTime
         T = title(sprintf('$t=%.2f$\n Arrived at time $t^* \\approx %.2f$',(n-1)*dt,endTime)); T.FontSize = 16; T.Interpreter = 'latex';
    end
    if PRINT
        % print images
        NAME = sprintf('%s%i',PRINTNAME,n-1);
        print(NAME,'-dpng');
    else
        pause(0.2);
    end
end

%%%%% FOR LAST PIC FROM PAPER
% print('ex2_pic5','-dpng');
%%%%%

end

