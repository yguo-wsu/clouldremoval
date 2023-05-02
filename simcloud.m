function C = simcloud(rows,cols) 
% This function simulate the cloud image of size rows and cols 
% rows: number of rows of the simulated image (>>64)
% cols: number of columns of the simulated image (>>64)
% It uses Perlin noise and start from a 64x64 simulation and then resize
% to the specified size due to the slowness of Perlin noise generator, 
% % which is the best generator. 

n = 64;
m = 64;
C = zeros(n, m);
C = perlin_noise(C);

% figure; imagesc(im); colormap gray;
C = imresize(C,[rows cols]);
C = imadjust(C,[],[],4); 
C = adapthisteq(C);
C(C<0) = 0;
C(C>1) = 1; 
end 
