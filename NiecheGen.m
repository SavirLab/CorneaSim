function [ Stem_cell  ] = NiecheGen( r,n,grid_circular,inter_ind )
%Get n number of cells postions inter_ind
%   Generate for ech stem cell niche vector includind 2 adjacent stem cells
%   and evreything in betwwen using the circular grid
grid = zeros(2*r+3,2*r+3);
[I,J]=ind2sub(size(grid),inter_ind);
Pos_vec=[I;J];
Pos_vec=Pos_vec';

D = pdist(Pos_vec);
S =squareform(D);

for i =1:n% S is nXn
    S_line=S(i,:);
    [~,IX] = sort(S_line);
    Stem_cell{1,i}=sub2ind(size(grid),Pos_vec(IX(1),1),Pos_vec(IX(1),2));
    Stem_cell{3,i}=sub2ind(size(grid),Pos_vec(IX(2:3),1),Pos_vec(IX(2:3),2));
end

%add niches to cells
C_cells=find(grid_circular==1);
[x_C,y_C]=ind2sub(size(grid),C_cells);
for i=1:length(C_cells)
    D=sqrt((Pos_vec(:,1)-x_C(i)).^2+(Pos_vec(:,2)-y_C(i)).^2);
    [~,IX]=sort(D);
    S_fix=inter_ind(IX(1:2));
    for k=1:2
        S_C_ind=find([Stem_cell{1,:}]== S_fix(k));
        Stem_cell{2,S_C_ind}=unique([Stem_cell{2,S_C_ind},C_cells(i)]);
    end
end



%remove unwanted cells
grid3=zeros(2*r+3,2*r+3);
for t=1:length(Stem_cell)
    grid3(Stem_cell{2,t})= grid3(Stem_cell{2,t})+1;
end
%imagesc(grid3)
prob_cells=find(grid3==1)';
[x_prob,y_prob]=ind2sub(size(grid),prob_cells);
for i=1:length(prob_cells)
    D=sqrt((Pos_vec(:,1)-x_prob(i)).^2+(Pos_vec(:,2)-y_prob(i)).^2);
    [~,IX]=sort(D);
    S_fix=inter_ind(IX(1:2));
    for k=1:2
        S_fix_ind=find([Stem_cell{1,:}]== S_fix(k));
        Stem_cell{2,S_fix_ind}=unique([Stem_cell{2,S_fix_ind},prob_cells(i)]);
    end
end

for i=1:n
    Stem_cell{3,i}=Stem_cell{3,i}';
    Stem_cell{2,i}=intersect(find(grid_circular==1)',Stem_cell{2,i});
    Stem_cell{2,i}=Stem_cell{2,i}(Stem_cell{2,i}~=Stem_cell{1,i});
    Stem_cell{2,i}=Stem_cell{2,i}(Stem_cell{2,i}~=Stem_cell{3,i}(1));
    Stem_cell{2,i}=Stem_cell{2,i}(Stem_cell{2,i}~=Stem_cell{3,i}(2));
end


%draw result 
%{
figure(2);hold on
grid2 = zeros(2*r+3,2*r+3);
color2=grid2;
color2(inter_ind)=100;
pause
for k=1:length(Stem_cell(1,:))
    color2(Stem_cell{2,k})=150;
    imagesc(color2)
    pause
    color2(Stem_cell{3,k})=200;
    color2(Stem_cell{1,k})=230;
    imagesc(color2)
    pause
    color2=grid2;
    color2(inter_ind)=100;
end 
imagesc(color2)

color3=color2;
 for S=1:length(Stem_cell)
    color2([Stem_cell{2,S}])=310;
    imagesc(color2)
    pause
    color2=color3;
end
%}
