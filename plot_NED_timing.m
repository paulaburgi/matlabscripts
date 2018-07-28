close all

bbox_ll=[-125 41;-120 46];
a=shaperead('FESM_13','BoundingBox',bbox_ll);

figure
for i=1:length(a)
    patch([a(i).X(1:end-1)],[a(i).Y(1:end-1)],a(i).s_date)  % the nan at the end of each x,y vector messes up the patch command.
end
axis image
caxis([1950 2017])
colormap jet
colorbar 



hold on; 
f = [-123.438 43.909; -123.299 43.909; -123.299 43.769; -123.438 43.769; -123.438 43.909];
plot(f(:,1), f(:,2), 'k', 'linewidth', 2); 

