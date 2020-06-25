 function plot_grid(s);

[ny,nx]=size(s);

hold on;
for ix = 1:nx; 
for iy = 1:ny; 
if s(iy,ix)==1; plot(ix,iy,'.');  end
end;
end;
