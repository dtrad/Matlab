I = imread('Water-Lilies.jpg');
 class_of_I = class(I);
 [x y] = meshgrid(1:256);
 [xi yi] = meshgrid(1:0.1:256);
 
 New_Image = cast(interp2(x,y,double(I),xi,yi,'linear'),class_of_I);