%% Define bouderies of rectangle lattice
function [ind,ind2,ind3,grid] = CircleGrid(r)

r2 = r-1;

theta = [-pi:0.00001:pi];
x = round(r*cos(theta))+r+2;
y = round(r*sin(theta))+r+2;

x2 = round(r2*cos(theta))+r+2;
y2 = round(r2*sin(theta))+r+2;



%figure;hold on
%plot(x,y)
%plot(x2,y2,'r')

grid = zeros(2*r+3,2*r+3);
ind = unique(sub2ind(size(grid), x, y));
ind2  = unique(sub2ind(size(grid), x2, y2));

grid(ind) = 1;
grid(ind2) = 2;

ind= find(grid==1);
ind=ind';

bw = im2bw(~grid);
CC = bwconncomp(bw);
for i = 1:length(CC.PixelIdxList)
    ind_cc{i} = CC.PixelIdxList{i};
    size_cc(i) = length(CC.PixelIdxList{i});
end
[val,ind_max] = max(size_cc);
ind3 = CC.PixelIdxList{ind_max};


grid(ind3) = 3;

% figure
% imagesc(grid);