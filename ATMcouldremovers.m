function res = ATMcouldremovers(imData,method,lambda,ratio,verbose)
% Cloud removal by ATM/RPCA on single channel images sequence
% imData: input image data, a 3D array with 3rd mode being the time
% method: which method to use (1-4). 
% lambda: regularisation parameter. Can be set by the equation used in the
% paper. 
% ratio: a tuning parameter to generate cloud mask. Default 0.8.
% verbose: 0/1 for printing out. Default 0. 

[rows,cols,n] = size(imData);

if ~exist('method','var')
    method = 1;
end

if ~exist('lambda','var')
    lambda = 1/sqrt(rows*cols);
end

if ~exist('ratio','var')
    ratio = 0.8; % this can be tuned
end

if ~exist('verbose','var')
    verbose = 0;
end

res.p1result=zeros(rows,cols,n);
res.clouds = res.p1result;
res.mask = false(rows,cols,n);
D=reshape(imData,rows*cols,n);

if verbose
    disp('Cloud detection')
end

tic 
switch method
    case 1
        res.method = 'Nonnegative RPCA';
        [B,E]=inexact_alm_r1pca_YG(D,lambda);
        %         [B,E]=inexact_alm_rpca(D,lambda);
    case 2
        res.method = 'ATM';
        [B,E]=bcs_exact(D,lambda,1e-5,100);
    case 3 
        res.method = 'ATM(RPCA+small noise)';
        [B,E,N]=rpca4atm(D,lambda);
        res.noise = reshape(N,rows,cols,n);
    case 4
        res.method = 'RPCA';
        [B,E]=inexact_alm_rpca(D,lambda);
end
res.optimisationtimeused = toc; 
res.p1result=reshape(B,rows,cols,n);
res.clouds=reshape(E,rows,cols,n);
for j=1:n
    res.mask(:,:,j)=abs(res.clouds(:,:,j))>ratio*std(E(:));
end

%morphological filters
se=strel('disk',3);
res.mask=imopen(res.mask,se);
res.mask=imdilate(imdilate(res.mask,se),se);

%cloud and shadow removal by matrix completion
[B_,E_]=inexact_alm_rpca_BS(D,reshape(res.mask,rows*cols,n));
res.result=reshape(B_,rows,cols,n);
res.totaltimeused = toc; 
