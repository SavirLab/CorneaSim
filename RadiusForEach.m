function [ rad_cell,rad_cell_half,rad_cell_best_final ] = RadiusForEach( r,grid,r_around )
%for each cell make availible indexes
% [ rad_cell,rad_cell_best,rad_cell_best2,rad_cell_best_n,rad_cell_best_n2, rad_cell_best_final ] = RadiusForEach( r,grid,r_around )
rad_cell={};
rad_cell_best={};
[rows,cols]=size(grid);

n_best=r_around*3;%best of 5 closest to the center
theta = [-pi:0.00001:pi];
x=[];
y=[];
for i=0:0.3:r_around
    x=[x,round(i*cos(theta))+r];
    y =[y,round(i*sin(theta))+r];
end
linearInd2= unique(sub2ind([rows cols], x, y));
[I,J] = ind2sub([rows cols],linearInd2);
V2=[I'-r,J'-r];
V2(find(~any(V2')),:)=[];
grid_temp=find(grid);

for i=1:rows
    for j=1:cols
        V_i=[V2(:,1)+i,V2(:,2)+j];
        V_i(:,1)=max( V_i(:,1),1);
        V_i(:,1)=min( V_i(:,1),rows);
        V_i(:,2)=max( V_i(:,2),1);
        V_i(:,2)=min( V_i(:,2),rows);
        V_check=sqrt((V_i(:,1)-r).^2 +(V_i(:,2)-r).^2);
        
        [~,sortIndex] = sort(V_check(:),'ascend');  %# Sort the values in
        maxIndex = sortIndex(1:n_best);  %# Get a linear index into A of the 5 largest values
        maxIndex_half = sortIndex(1:round(length(V_check)/2));  %# get half circle closes to center
        V_half=V_i(maxIndex_half,:);
        V_to_center=[r-i,r-j];
        V_half_center=r-V_half;
        V_to_center=V_to_center/norm(V_to_center);
        for iii=1:length(V_half)
        norm_V_half(iii)=norm(V_half_center(iii,:));
        end 
        costheta= V_to_center*V_half_center'./norm_V_half;
        
        [~,ind_c] = max(costheta);
        [~,ind_c_2] = min(V_check);
        [~,sorttheta] = sort(costheta(:),'descend');
        
        sorttheta_ind=sorttheta(1:n_best);
        
        sorttheta_ind_half=sorttheta(1:round(length(sorttheta)/2.8));
        V_i_best=V_half(ind_c,:);
        V_i_best2=V_i(ind_c_2,:);
        V_i_best_n=V_half(sorttheta_ind_half,:);
        V_i_best_n2=V_i(maxIndex,:);
        
        V_ind=sub2ind([rows cols], V_i(:,1), V_i(:,2));
        V_ind_half=sub2ind([rows cols], V_half(:,1), V_half(:,2));
        %V_ind_best=sub2ind([rows cols], V_i_best(:,1), V_i_best(:,2));
        %V_ind_best2=sub2ind([rows cols], V_i_best2(:,1), V_i_best2(:,2));
        V_ind_best_n=sub2ind([rows cols], V_i_best_n(:,1), V_i_best_n(:,2));
        V_ind_best_n2=sub2ind([rows cols], V_i_best_n2(:,1), V_i_best_n2(:,2));
        rad_cell{i,j} = intersect(V_ind,grid_temp);
        rad_cell_half{i,j}=intersect(V_ind_half,grid_temp);
        %rad_cell_best2{i,j}=intersect(V_ind_best2,grid_temp);
        rad_cell_best_n{i,j}=intersect(V_ind_best_n,grid_temp);
        rad_cell_best_n2{i,j}=intersect(V_ind_best_n2,grid_temp);
        rad_cell_best_final{i,j}=intersect(rad_cell_best_n{i,j},rad_cell_best_n2{i,j});
    end
end



%squre ind
%{
vxn = r_around:-1:-r_around;
[VXN,VXY] = meshgrid(vxn,vxn);
V = [VXN(:),VXY(:)];
V(find(~any(V')),:)=[];
B=V+50;
linearInd = sub2ind([rows cols], B(:,1), B(:,2));
%}

end
%%check
%{
[rows,cols]=size(grid);
color=grid;
i=60;
j=90;
for k=1:40
    color(r,r)=40;
    color(i,j)=10;
    color(rad_cell_best{i,j})=20;
    color(rad_cell_best_final{i,j})=30;
    %color(rad_cell_best_n2{i,j})=40;
    %color(rad_cell_best{i,j})=30;
    %color(rad_cell_best2{i,j})=40;
    imagesc(color)
    G=randi([1,length(rad_cell_best_final{i,j})],1);
    [i,j]=ind2sub([rows,cols],rad_cell_best_final{i,j}(G));
end
%}