function [ inner_ind ] = GenInCirc(r,grid_boundery,bias_area)
%Generate circle with 1/3 of the area of the big circle in the center

grid=grid_boundery.*0;

r2 = sqrt(bias_area)*r;
theta = [-pi:0.00001:pi];

x = round(r2*cos(theta))+r+2;
y = round(r2*sin(theta))+r+2;
inner_ind = unique(sub2ind(size(grid), x, y));

grid(inner_ind)=1;
bw = im2bw(~grid);
CC = bwconncomp(bw);
for i = 1:length(CC.PixelIdxList)
    ind_cc{i} = CC.PixelIdxList{i};
    size_cc(i) = length(CC.PixelIdxList{i});
    is1=any(CC.PixelIdxList{i}==1);
    if is1
        ind_around=CC.PixelIdxList{i};
    end
end
[val,ind_max] = max(size_cc);
%ind3 = CC.PixelIdxList{ind_max};
ind3 = ind_around;
grid=ones(size(grid));

ind3=ind3';
grid(ind3)=0;
inner_ind=find(grid==1);

end

