function [ trajec_mat ] = Pushgrid(r,grid, x_new,y_new,x_hole,y_hole,cent,phi_cell)
%push cells


%% Test allign holes in rect gri
%
% [ind1,ind2,ind3,grid_boundery] = CircleGrid(r);% grid(1:L/Ns:end,1)=0;
% 
% grid = grid_boundery;
% grid(ind1) = 1;
% grid(ind2) = 1;
% grid(ind3) = 1;


%%create grid 
[~,cols]=size(grid);
%[XX,YY] = meshgrid([1:rows],[1:cols]);

%ind_1 = find(grid==1);

%grid_vec = grid(:);


%xx = XX(:);
%yy= YY(:);
  %}
%plot(xx,yy,'ko');hold on;
%plot(xx(grid==1),yy(grid==1),'bo')


%x_hole = xx(ind_d);
%y_hole = yy(ind_d);

%ind_d=sub2ind([rows,cols],x_hole,y_hole);
%xx(ind_d) = [];
%yy(ind_d) = [];

%theta = 2*pi*rand(size(xx));

% xx = xx+0.1*cos(theta);
% yy = yy+0.1*sin(theta);

%xx = [xx; xx_new];
%yy = [yy ;yy_new];

%figure(2)
%plot(xx,yy,'ko');hold on
%figure(2)
%plot(xx_new,yy_new,'ro');
%plot(x_hole,y_hole,'bo');

% find unit vector between the hole and the new born

%%allocate for speed calculate meshgrid and V_normal

vxn = [1 0 -1];
vyn = [1 0 -1];

[VXN,VXY] = meshgrid(vxn,vyn);
V = [VXN(:),VXY(:)];

%theta = 2*pi*rand(1);
%xx_new = x_new+0.5*cos(theta)';
%yy_new = y_new+0.5*sin(theta)';


v_c=[r r];
norm_v_c=norm(v_c);
if cent
    %old way chck who is closer
    %V_check=sqrt((-(V(:,1)+x_new)+r).^2 +(V(:,2)+y_new-r).^2);
    %[~,ind_c] = min(V_check);
    %xx_new = x_new+V(ind_c,1);
    %yy_new = y_new+V(ind_c,2);
    
    theta =  phi_cell(1)+phi_cell(5)*rand(1);
    theta(theta<-pi)=theta(theta<-pi)+2*pi;
    theta(theta>pi)=theta(theta>pi)-2*pi;
    xx_new = x_new+0.5*cos(theta)';
    yy_new = y_new+0.5*sin(theta)';
else
    theta = 2*pi*rand(1);
    xx_new = x_new+0.5*cos(theta)';
    yy_new = y_new+0.5*sin(theta)';
end


for iii=1:length(V)
    norm_V(iii)=norm(V(iii,:));
end 
V_int = V;
for k = 1:5*r
 vx = (xx_new-x_hole);
 vy = (yy_new-y_hole);
 v = [vx,vy];

 % vectors to nearest niebers to hole
 
%vxn = [1 0 -1];
%vyn = [1 0 -1];

%[VXN,VXY] = meshgrid(vxn,vyn);

V = V_int;
norm_V_run=norm_V;
V_del=[V(:,1)+x_hole,V(:,2)+y_hole];%find all moves
del_ind=[];

for inddd=1:length(V)%find  illigal moves
    V_d=V_del(inddd,:);
    V_d(1)=max(1,V_d(1));
    V_d(2)=max(1,V_d(2));
    V_d(1)=min(cols,V_d(1));
    V_d(2)=min(cols,V_d(2));
    if grid(V_d(1),V_d(2))== 0%if its zero in grid mat its not legal 
        del_ind=[del_ind,inddd];
    end
end

V(del_ind,:)=[];%delete illigal moves
norm_V_run(del_ind)=[];


t_s=(v*V')./(norm(v).*norm_V_run);%matrix code for vector dot product
theta2 = acos(t_s);
[~,ind] = min(theta2);

x_hole_old=x_hole;
y_hole_old=y_hole;

%xx = [xx; x_hole];
%yy = [yy; y_hole];

trajec_mat{1,k}=[x_hole,y_hole];
trajec_mat{2,k}=[x_hole + V(ind,1),y_hole + V(ind,2)];
x_hole = x_hole + V(ind,1);
y_hole = y_hole + V(ind,2);

%tempx = (xx==x_hole);
%tempy = (yy==y_hole);
%tempxy = find(tempx.*tempy);
%tempxy=find(xx==x_hole & yy==y_hole);


%xx(tempxy) = [x_hole_old];
%yy(tempxy) = [y_hole_old];
%theta=[];
%t=[];
t_s=[];
theta2=[];

%  figure(2);hold on
%  plot(xx,yy,'ko');hold on
%  figure(2);hold on
%  plot(xx_new,yy_new,'ro');

% figure(2)
% plot(x_hole,y_hole,'bo');

% stop criteria

%1 if the distance between the new hole and the new born is smaller than 1

if sqrt( (xx_new-x_hole)^2+(yy_new-y_hole)^2 )<1
    
   %xx = [xx; x_hole];
    %yy = [yy; y_hole];
    
    %figure(k+3)
    %plot(xx,yy,'ko');hold on
    %figure(k+3)
    %plot(xx_new,yy_new,'ro');
    if 0
        figure(666)
        [~,t_cols]=size(trajec_mat);
        for ttt=1:t_cols
            color(trajec_mat{2,ttt}(1),trajec_mat{2,ttt}(2))=400;
        end
        color(trajec_mat{1,1}(1),trajec_mat{1,1}(2))=700;
        color(i,j)=600;
        figure(66)
        imagesc(color)
    end
    break
end
%figure(2)
end
% % sort according to y
% 
% [yy_sort,ind_yy_sort] = sort(yy);
end
