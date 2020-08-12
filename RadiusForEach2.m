function [ rad_cell ] = RadiusForEach2( r,grid,r_around,phi_cell)
%for each cell make availible indexes
% [ rad_cell,rad_cell_best,rad_cell_best2,rad_cell_best_n,rad_cell_best_n2, rad_cell_best_final ] = RadiusForEach( r,grid,r_around )


%[ phi_cell] = PhiForEach( r,grid,0.2*pi );% Debug1

rad_cell={};% for each cell list of possibile outcomes
[rows,cols]=size(grid);

theta = [-pi:0.00001:pi];
x=[];
y=[];
for i=0:0.3:r_around
    x=  [x,round(i*cos(theta))+r];
    y =[y,round(i*sin(theta))+r];
end
linearInd2= unique(sub2ind([rows cols], x, y));
[I,J] = ind2sub([rows cols],linearInd2);%find all possible combinations 
V2=[I'-r,J'-r];
V2(find(~any(V2')),:)=[];%remove [0 0]
grid_temp=find(grid);%circle   indexes 

[V_phi,~]= cart2pol(V2(:,1),-1.*V2(:,2));% find all the direction angels 

for i=1:rows
    for j=1:cols
        V_i=[V2(:,1)+i,V2(:,2)+j];%V_i=[x,y]+[all possible angles]
        %delte non grid results(edges)
        V_i(:,1)=max( V_i(:,1),1);
        V_i(:,1)=min( V_i(:,1),rows);
        V_i(:,2)=max( V_i(:,2),1);
        V_i(:,2)=min( V_i(:,2),rows);
        
        phi_MAT=phi_cell{i,j};%phi cell =[theta_min,theta_max,phi,bigopening,dPhi];
        %phi_MAT3(i,j)=phi_MAT(3);
        
        DEL_phi_ind=unique([find(V_phi>phi_MAT(2))',find(V_phi<phi_MAT(1))']);%find all angles bigger than max phi OR smaller than min phi
        
        phi_del=V_phi;
        phi_del(DEL_phi_ind)=[];
        V_i_del1=V_i(DEL_phi_ind,:);
        V_i_del2=V_i;
        V_i_del2(DEL_phi_ind,:)=[];%delete angles with DEL_phi_ind = V_i_del2 result ! 
        
        %{
        if ((phi_MAT(3)-max(phi_del)<0.0476) && (phi_MAT(3)-min(phi_del))>=-0.0476)...
                || (max(phi_del)==min(phi_del))
                
            v_final=V_i_del2;
        else
            v_final=V_i_del1;
        end
        %}
        if length(phi_del)<(length(V2)/2) %choose always the small portion of the circle
            v_final=V_i_del2;
        else
            v_final=V_i_del1;
        end

        if phi_MAT(4)%if dPhi >= 3.14 choose the portion with the angle inside it . 
% 
%             if length(phi_del)>(length(V2)/2) 
%             v_final=V_i_del2;
%             else
%             v_final=V_i_del1;
%             end

            if isempty(V_i_del2)
                v_final=V_i_del1;
            elseif ((phi_MAT(3)-max(phi_del)<0.0476) && (phi_MAT(3)-min(phi_del))>=-0.0476)...
                    %|| (max(phi_del)==min(phi_del))

                v_final=V_i_del2;
            else
                v_final=V_i_del1;
            end
        end
        
        V_ind=sub2ind([rows cols], v_final(:,1), v_final(:,2));

        rad_cell{i,j} = intersect(V_ind,grid_temp);
        

    end
    
    

end



        %debug
%{
    DB=cellfun('length',rad_cell)
imagesc(DB)

[rows,cols]=size(grid);
color=grid;

for i=195
    for j=75
%i=60;%140 100!!!
%j=140;
    for k=1:1
    color(r+2,r+2)=40;
    color(i,j)=10;
    color(rad_cell{i,j})=20;
    %color(rad_cell_best_final{i,j})=30;
    %color(rad_cell_best_n2{i,j})=40;
    %color(rad_cell_best{i,j})=30;
    %color(rad_cell_best2{i,j})=40;
    
    if length(rad_cell{i,j})>1
        G=randi([1,length(rad_cell{i,j})],1);
        [i,j]=ind2sub([rows,cols],rad_cell{i,j}(G));
    elseif length(rad_cell{i,j})==1
        [i,j]=ind2sub([rows,cols],rad_cell{i,j});
    end
    
    end
    end
end
imagesc(color)
print('ok')
        %}
