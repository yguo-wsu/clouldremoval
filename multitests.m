%%
% nsims = 2; nrepeats=2; imgpath = './tyrol-w1.tif'; 
% Use the following linux script to run:
% nohup matlab -nosplash -noFigureWindows -r "try; nsims=7; nrepeats=50;imgpath='./kitsap1.tif'; multitestsontaichi; catch; end; quit" > out7 &

% Set default values
if ~exist('imgpath','var')
    imgpath='./tyrol-w1_1k.tif';
end

if ~exist('nsims','var')
    nsims=3;
end

if ~exist('nrepeats','var')
    nrepeats = 2;
end

[~, compname] = system('hostname');
compname = compname(find(~isspace(compname))); % This is a must. There are white spaces everywhere. 

[~, imgfn]=fileparts(imgpath);

% fn = ['tmpmultitests_',imgfn,'_nsim=',num2str(nsims),'_',char(datetime('now','TimeZone','local','Format','yMMddHHmmss')),'.mat'];

im = imread(imgpath);
maxDN = 2^ceil(log2(double(max(im(:))))); 
disp(['The image max DN is ',num2str(maxDN)])

im = double(imresize(rgb2gray(im), [1024,1024]))/maxDN;
[rows, cols] = size(im);

disp(['The image is of size ',num2str(rows), ' by ', num2str(cols)])

lambda = 1/sqrt(rows*cols);
lambdarange = exp(linspace(log(lambda/10),log(lambda*10),51));
lambdarange(26) = lambda;
lambdarange = lambdarange(5:27);
nparas = length(lambdarange);
ntotal = nrepeats*nparas;

R = zeros(nrepeats*10,nparas);
time_total = zeros(nrepeats*4,4);
time_optimisation = time_total;

jobname = ['multitests_',imgfn,'_nsim',num2str(nsims),'_nrepeats',num2str(nrepeats)];
disp(['I am running the job ',jobname])

t1=tic;
parfor ipara=1:nparas
    R_ = zeros(nrepeats,10);
    time_total_ = zeros(nrepeats,4);
    time_optimisation_ = zeros(nrepeats,4);
    t2=tic;
    for irep=1:nrepeats
        k = (ipara-1)*nrepeats+irep;
        imData = zeros(rows, cols, nsims);
        cData = zeros(rows, cols, nsims);
        for i=1:nsims
            c = simcloud(rows,cols);
            imData(:,:,i) = c + (1-c).*im;
            cData(:,:,i) = c;
        end
        lambda = lambdarange(ipara);
        rpcann = ATMcouldremovers(imData,1,lambda);
        atm = ATMcouldremovers(imData,2,lambda);
        atm2 = ATMcouldremovers(imData,3,lambda);
        rpca = ATMcouldremovers(imData,4,lambda);
        R_(irep,:) = [irep, lambda, evalfunc(atm2.result,im), ...
            evalfunc(atm.result,im), evalfunc(rpcann.result,im), ...
            evalfunc(rpca.result,im), evalfunc(atm2.p1result,im), ...
            evalfunc(atm.p1result,im), evalfunc(rpcann.p1result,im), ...
            evalfunc(rpca.p1result,im)];
%         time_optimisation(k,:) = [atm2.optimisationtimeused, atm.optimisationtimeused, rpcann.optimisationtimeused, rpca.optimisationtimeused];
%         time_total(k,:) = [atm2.totaltimeused, atm.totaltimeused, rpcann.totaltimeused, rpca.totaltimeused];
        time_optimisation_(irep,:) = [atm2.optimisationtimeused, atm.optimisationtimeused, rpcann.optimisationtimeused, rpca.optimisationtimeused];
        time_total_(irep,:) = [atm2.totaltimeused, atm.totaltimeused, rpcann.totaltimeused, rpca.totaltimeused];
        fprintf("lambda=%f, trial %d done. Total: %f\n",lambda,irep, k/ntotal);
    end
    R(:,ipara) = R_(:);
    time_total(:,ipara) = time_total_(:);
    time_optimisation(:,ipara) = time_optimisation_(:);

%     eval(['save ',fn, ' R time_optimisation time_total'])
    msg=[compname, ' tested the ', num2str(ipara), 'th parameter out of ', num2str(nparas), ' for ', jobname, ' for ', num2str(toc(t2)), 'seconds.'];
    disp(msg)
%     msg=['echo "',msg,'"| mail -s ', compname, ' Progress report y.guo@uws.edu.au'];
%     system(msg);
end

R = reshape(permute(reshape(R,nrepeats,10,nparas),[1,3,2]),ntotal,10);
time_optimisation = reshape(permute(reshape(time_optimisation,nrepeats,4,nparas),[1,3,2]),ntotal,4);
time_total = reshape(permute(reshape(time_total,nrepeats,4,nparas),[1,3,2]),ntotal,4);


multisim.Rcolnames = {'indrepeat','lambda','ATM2+MC','ATM+MC','NNRPCA+MC','RPCA+MC','ATM2 ','ATM','NNRPCA','RPCA'};

fn = ['multitests_',compname,'_',imgfn,'_nsim=',num2str(nsims),'_',char(datetime('now','TimeZone','local','Format','yMMddHHmmss')),'.mat'];
multisim.R = R;
multisim.time_optimisation = time_optimisation;
multisim.time_total = time_total;
multisim.lambdas = lambdarange;
multisim.nsims = nsims;
multisim.nrepeats = nrepeats;
multisim.imgpath = imgpath;
multisim.compname = compname;
multisim.grosscomputationtime = toc(t1);
eval(['save ',fn, ' multisim'])

msg = ['Job ', jobname, ' has done!'];
disp(msg)
disp(['I used ',num2str(multisim.grosscomputationtime), ' seconds in total for this job.']);
msg = ['echo "',msg,'" | mail -s ',compname, ' for ', jobname, ' cloud removal tests done! y.guo@uws.edu.au'];
system(msg);

if ~ismac
    quit
end

%% Visualise results
multisim.Rcolnames = {'indrepeat','lambda','aATM+MC','ATM+MC','NNRPCA+MC','RPCA+MC','aATM ','ATM','NNRPCA','RPCA'};

[~,ind] = sort(multisim.R(:,2));
R = multisim.R(ind,:);
lambdas = unique(R(:,2));

visres(R(:,3:end),multisim.nrepeats,lambdas);
legend(multisim.Rcolnames(3:end))

xline(1/2^10)

%%
visres(multisim.time_optimisation,multisim.nrepeats,multisim.lambdas);
legend(multisim.Rcolnames(3:3+4-1))
xline(1/2^10)


%%
function visres(r,n,x,shift)
addpath('/Users/yiguo/FromownCloud/Yi/research projects/tools/distinguishable_colors/')
N = size(r,2);
cols = distinguishable_colors(N);
markers = ['o','+','*','.','x','_','|','s','d','^','v','>','<','p','h'];
figure, hold on
for i=1:N
    [m{i} s{i}] = getstat(r(:,i),n);
    %     plot(m{i}, 'Color',cols(i,:),'Marker',markers(i));
    errorbar(x,m{i},s{i}, 'Color',cols(i,:),'Marker',markers(i));
end
set(gca,'xscale','log')
end

function [m,s] = getstat(x,n)
x = reshape(x,n,[]);
m = mean(x,1);
s = std(x);
end

function x = stackimgseq(a,b)
x = a;
x(:,:,end+1:end+size(b,3)) = b;
end

function x=farray(a,n)
if size(a,3)==1
    x = repmat(a,[1,1,n]);
else
    x=a;
end
end

function v = evalfunc(result,im)
nsims = size(result,3);
v = norm(result(:)-repmat(im(:),nsims,1),'fro')/norm(repmat(im(:),nsims,1),'fro');
end