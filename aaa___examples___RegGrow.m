function aaa___examples___RegGrow
%run a few examples and tests

%% generate an example image
if exist('OCTAVE_VERSION', 'builtin')
    IM = get(0,'DefaultImageCData');IM = IM/max(IM(:));
else
    defimage = pow2(get(0,'DefaultImageCData'),47);
    IM = bitshift(defimage,-37);IM = fix(IM);
    IM = bitand(IM,31);IM = IM/max(IM(:));
end

%% show default
mkfig(1,'show default')
subplot(1,2,1),imshow(IM),title('original example image')
subplot(1,2,2),imshow(fuse_mask(IM,RegGrow)),title('segmented background')

%% showcase leak detection
%this only works with the Matlab example image
mkfig(2,'eye segmentation with and without leak detection')
opts=struct;
opts.seed=[24 30];
opts.maxDiff=0.15;

repeats=30;%repeat timing to account for the JIT and random variation
tic
for repeats=1:repeats
    r1=RegGrow(IM,opts);
end
t1=toc;

opts.overflow=100;
tic
for repeats=1:repeats
    r2=RegGrow(IM,opts);
end
t2=toc;

opts.overflow=[];opts.maxDiff=0.1;
tic
for repeats=1:repeats
    r3=RegGrow(IM,opts);
end
t3=toc;

IM1=fuse_mask(IM,r1);
IM2=fuse_mask(IM,r2);
IM3=fuse_mask(IM,r3);

subplot(2,2,1),imshow(IM),title({'original example image','( * = seed)'})
hold on,plot(opts.seed(2),opts.seed(1),'y*'),hold off
subplot(2,2,2),imshow(IM1)
title({'leaked segmentation',sprintf('(t=%.1fms)',1000*t1/repeats)})
subplot(2,2,3),imshow(IM2)
title({'with leak detection',sprintf('(t=%.1fms)',1000*t2/repeats)})
subplot(2,2,4),imshow(IM3)
title({sprintf('better maxDiff: %.2f',opts.maxDiff),...
    sprintf('(t=%.1fms)',1000*t3/repeats)})

end
function mkfig(num,titletxt)
f=figure(num);
clf(num)
set(f,'NumberTitle','off')
set(f,'Name',titletxt)
end
function fused=fuse_mask(IM,mask)
r_factor=0.5;
[R,G,B]=deal(IM,IM,IM);
R(mask)=r_factor+R(mask)*(1-r_factor);
G(mask)=G(mask)/2;
B(mask)=B(mask)/2;
fused=cat(3,R,G,B);
end