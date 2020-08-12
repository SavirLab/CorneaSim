function [ phi_cell] = PhiForEach( r,grid,dPhi )
%for each cell make availible
%[theta_min,theta_max,phi,bigopening,dPhi];
dPhi=dPhi/2;
phi_cell={};
[rows,cols]=size(grid);
for i=1:rows
    for j=1:cols
        x=i;
        y=cols-j;
        V=[r+2-x,r+2-y];
        [phi,~]=cart2pol(V(1),V(2));
        theta_min=phi-dPhi;
        theta_min(theta_min<-pi)=theta_min(theta_min<-pi)+2*pi;
        theta_max=phi+dPhi;
        theta_max(theta_max>pi)=theta_max(theta_max>pi)-2*pi;
        bigopening=dPhi>=pi/2;
        
        if theta_max<theta_min
            [theta_max,theta_min]=deal(theta_min,theta_max);
        end
        phi_cell{i,j}=[theta_min,theta_max,phi,bigopening,dPhi*2];
        
    end
end
