function [h] = showEmbryo(V,X,Y,Z,R,n,t,a)
    V(V<t)=NaN;
    h = slice(X,Y,Z,V,R.XWorldLimits(1):1:R.XWorldLimits(2),R.YWorldLimits(1):1:R.YWorldLimits(2),R.ZWorldLimits(1):1:R.ZWorldLimits(2));
    set(h, 'FaceColor', 'interp', 'EdgeColor', 'none','FaceAlpha',a, 'EdgeAlpha',a);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis equal;
    axis vis3d;
end