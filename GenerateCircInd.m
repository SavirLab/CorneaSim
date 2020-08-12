function [ind_cell,inter_ind,Stem_cell] = GenerateCircInd(r,n)
% generate index cell slicing radius r circle into n pieces
%Stem_cell= for each niche list of stem cells
r2=r-1;
[ind1,ind2,ind3,grid_boundery] = CircleGrid(r);
grid_circular = grid_boundery;
grid_circular(ind1) = 1;
grid_circular(ind2) = 0;
grid_circular(ind3) = 0;


theta = [-pi:0.00001:pi];
grid = zeros(2*r+3,2*r+3);
len_thet=length(theta);
b=0;
jump=round(len_thet/n);

rad_cell=RadiusForEach( r,grid_circular,1);

for i=1:n
    if i==1
        x = round(r*cos(theta(1:jump)))+r+2;
        y = round(r*sin(theta(1:jump)))+r+2;
        x2 = round(r2*cos(theta(1:jump)))+r+2;
        y2 = round(r2*sin(theta(1:jump)))+r+2;
    else
        x = round(r*cos(theta(jump*(i-1):min(jump*(i),len_thet))))+r+2;
        y = round(r*sin(theta(jump*(i-1):min(jump*(i),len_thet))))+r+2;
        x2 = round(r2*cos(theta(jump*(i-1):min(jump*(i),len_thet))))+r+2;
        y2 = round(r2*sin(theta(jump*(i-1):min(jump*(i),len_thet))))+r+2;
    end
    grid = zeros(2*r+3,2*r+3);
    ind = unique(sub2ind(size(grid), x, y));
    ind2  = unique(sub2ind(size(grid), x2, y2));
    grid(ind) = 1;
    grid(ind2) = 2;

    ind= find(grid==1);
    ind=ind';
    ind=intersect(find(grid_circular==1)',ind);
    ind_cell{i} = ind;
    x_i = round(r*cos(theta(jump*(i-1)+1)))+r+2;
    y_i = round(r*sin(theta(jump*(i-1)+1)))+r+2;
    
    [x_ind, y_ind]=ind2sub(size(grid),ind);
    D=(sqrt((x_ind-x_i).^2+(y_ind-y_i).^2));
    %ind_stem=sub2ind(size(grid), x_i, y_i);
    [~,index] = min(D);
    inter_ind(i)=ind(index);
end




[ Stem_cell  ] = NiecheGen( r,n,grid_circular,inter_ind );
%old code
%{
% new find stem cell
figure(2);hold on
grid2 = zeros(2*r+3,2*r+3);
color2=grid2;
for k=1:length(ind_cell)
    color2(ind_cell{k})=k*10;
    imagesc(color2)
end
S_ind_prog=[]
for T=1:n
    for K=1:length(ind_cell{T})
        list=ind_cell{T};
        [i,j]  = ind2sub([size(grid)],list(K));
        %find diffrent colored cells
        cand=rad_cell{i,j}(color2(rad_cell{i,j})~=color2(list(K)));
        if ~isempty(cand)
            %S_ind(T)=cand;
            %break
            S_ind_prog=[S_ind_prog,cand']
        end
        %% TO DO
        1.Remove Sprog that are thoucing other Sprogs 
        2.for each intersection choose 1  sprog and save their S cell
        3. ??? ECH color 1 STEM CELL!!! 
        4.profit
        
    end
    
end




%Generate stem for each niche

Stem_cell={};
for S=1:length(inter_ind)
    Stem_cell{1,S}=[];
    Stem_cell{2,S}=[];
end
for S=1:length(inter_ind)
    for f=1:length(ind_cell) %search for the area to combine
        if ismember(inter_ind(S), ind_cell{f})
            Stem_cell{1,S}=inter_ind(S);
            Stem_cell{2,S}=[Stem_cell{2,S},ind_cell{f}];
            Stem_cell{2,S}=unique(Stem_cell{2,S});
        end
    end
    if length(intersect(Stem_cell{2,S},inter_ind))<3
        disp('now')
        if S<length(ind_cell)
            Stem_cell{2,S}=[Stem_cell{2,S},ind_cell{S+1}];
        else
            Stem_cell{2,S}=[Stem_cell{2,S},ind_cell{1}];
        end
        
    end
end
%}

X=cellfun(@(x) length(intersect(x,inter_ind)) ,Stem_cell(2,:));

if any(X>0)
    disp('ERROR');
end

%draw result 
%{
figure(2);hold on
grid2 = zeros(2*r+3,2*r+3);
color2=grid2;
color2(inter_ind)=100;
for k=1:length(ind_cell)
    color2(ind_cell{k})=k*10;
    imagesc(color2)
    %color(min(ind_cell{k}(1)))=100;
end
color2(S_ind_prog)=length(ind_cell)*10+50;
imagesc(color2)
color3=color2;
 for S=1:length(Stem_cell)
    color2([Stem_cell{2,S}])=310;
    imagesc(color2)
    pause
    color2=color3;
end
%}
