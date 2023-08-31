function [xC,r] = generateDisjointCircles(a,b,c,d,minRad,maxRad,DIST,howMany)
%
% This is a helper function for randomly generating disjoint circles you 
% can use as obstacles. Inputs:
%   - (a,b) x (c,d), the centers of the circles will be in this square
%   - minRad / maxRad, the minimum and maximum possible radius of each
%     circle
%   - DIST, the minimum distance between any two circles
%   - howMany, the number of circles you want
%
%   It will take 10000 random pulls to generate the circles. With certain
%   parameters, it will be impossible to find more circles; in that case,
%   it will generate as many disjoint circles as it can in 10000 tries.
%
xC = [];
r = [];
DIST = max(DIST,0);
k = 1;
while k <= 10000 && length(r)<howMany
   x = a + (b-a)*rand(1,1); y = c + (d-c)*rand(1,1);
   RAD = minRad + (maxRad - minRad)*rand(1,1);
   if isempty(r)
      r = RAD;
      xC = [x;y];
   else
       check = 1;
       for m = 1:length(r)
            check = check*(norm(xC(:,m) - [x;y],2)>=r(m)+RAD+DIST); 
       end
       if check
           xC(:,end+1) = [x;y]; r(end+1) = RAD;
       end
   end
   k=k+1;
end

if length(r) < howMany
   fprintf('Failed to generate %i disjoint circles with centers in [%.2f,%.2f] x [%.2f,%.2f] in %i trials.\n',howMany,a,b,c,d,k); 
else
   fprintf('Generated %i disjoint circles with centers in [%.2f,%.2f] x [%.2f,%.2f] in %i trials.\n',howMany,a,b,c,d,k);  
end
end

