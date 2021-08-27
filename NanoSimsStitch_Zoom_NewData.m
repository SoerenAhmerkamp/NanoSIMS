clear all
figure, 

% Load export from Look@NanoSIMS
files2 = dir('*.mat');

numb = [];
posx = [];
posy = [];
posz = [];
c12 = {};
n15 = {};
n15_int = {};
vlas = [18 19 20 21 22 23 77 78 79 80 81 82 83 84 88 89 90 91 92 93 94];
%vlas = [77 78 79 80 81 82 83 84 88 89 90 91 92];
%vlas = [1:18 20:25]; % Inner circle
o = 0;
for m = vlas
    o = o+1;
    load(files2(m).name,'p','im','planes')
    files2(m).name
    strT = strsplit(files2(m).name,'_');
    %strT = strsplit(strT{4},'.');
    strT = strsplit(strT{3},'.');
    
    numb(end+1) = str2num(strT{1});
    posy(end+1) = p.pos(1);
    posx(end+1) = p.pos(2);
    posz(end+1) = p.pos(3);
    
    for n = 1:length(im)
       im{n} = double(im{n});
    end
    
    imS{m} = im;
    c12{m} = sum(im{1},3);
    n15{m} = sum(im{4},3)./(sum(im{3},3));
    %n15{m} = sum(im{4},3)./(sum(im{3},3));
    %n15{m} = median(double(im{4})./double(im{3}),3);
    %n15(m) = median(n15_temp,3).*40;
    
    n15_int{m} = sum(im{3},3);
    
    %pause
    temp = im{1};
    %imagesc(temp(:,:,20))
end
%[numb,ind] = sort(numb);
%posx = posx(ind);
%posy = posy(ind);


%% Plot Positions

posx = posx-min(posx);
posy = posy-min(posy);

hold on
plot(posx,posy,'o');
for n = 1:length(posx)
    text(posx(n),posy(n),num2str(numb(n)))
end

dx = 256/20;

oposx = posx;
oposy = posy;
%%
% posx(14) = posx(14)+1;
% posx(13) = posx(13)+1;
% posx(12) = posx(12)+1;
% posx(11) = posx(11)+1;
% 
% posx(15) = posx(15)+1;
% posx(16) = posx(16)+1;
% posx(17) = posx(17)+1;
% 
% posx(21) = posx(21)+1;
% posx(20) = posx(20)+1;
% 
% posx(18) = posx(18)+1;
% posx(19) = posx(19)+1;
% 
% posx(21) = posx(21)+1;
% posx(20) = posx(20)+1;
% 
% posx(1) = posx(1)+1;
% posx(2) = posx(2)+1;
% posx(3) = posx(3)+1;
%%
posx = oposx;
posy = oposy;

stitchtMat = zeros(round( (max(posy)+20)*dx ),round((max(posx)+20)*dx) );
stitchtMat_int = zeros(round( (max(posy)+20)*dx ),round((max(posx)+20)*dx) );
tempMat = zeros(round( (max(posy)+20)*dx ),round((max(posx)+20)*dx),size(im{1},3));

imNew = {};
for j = 1:length(im)
   imNew{j} = zeros (round( (max(posy)+20)*dx ),round((max(posx)+20)*dx),size(im{1},3));
end

% Uncomment below for correction determined through crosscorr
%posx(5:8) = posx(5:8)+1.5;
%posy(1) = posy(1)-2;
%%%posy(1) = posy(1)-2;
%%%posy(4) = posy(4)-2;
%posx(9:11)=posx(9:11)+2;
%posy(9:11)=posy(9:11)+1;
%posx(10)=posx(10)-0.9;
%posy(11)=posy(11)-0.3;
%posx(11)=posx(11)-0.3;
%posy(5:8) = posy(5:8)+0.8;
%posx(12:13) = posx(12:13)+2.4;
%posy(12:13) = posy(12:13)+1.7;
%posx(13) = posx(13)-0.8;
%posy(13) = posy(13)+0.4;

for m = 1:length(posy)
    xInd = posy(m)*dx+1;
    yInd = posx(m)*dx+1;
    n = vlas(m);
    
    % if size(c12{n})
    stitchtMat(xInd+1:xInd+255-1,yInd+1:yInd+255-1) = n15{n}(2:end-1,2:end-1);
    stitchtMat_int(xInd+1:xInd+255-1,yInd+1:yInd+255-1) = n15_int{n}(2:end-1,2:end-1);
    
    imagesc(stitchtMat);
    caxis([0 0.04])
    
    for j = 1:length(im)
        for i = 1:size(im{j},3)
            imNew{j}(xInd+1:xInd+255-1,yInd+1:yInd+255-1,i) = imS{n}{j}(2:end-1,2:end-1,i);
        end
        %imagesc(imNew{j}(:,:,20))

    end
    
    pause(0.01);
    %end
end

%
v = stitchtMat_int;
colormap(gray)
map = colormap;
minv = min(v(:));
maxv = max(v(:));
ncol = size(map,1);
s = round(1+(ncol-1)*(v-minv)/(maxv-minv));
rgb_image = ind2rgb(s*5,map);
imshow(rgb_image)
%
%%
imagesc(stitchtMat, 'AlphaData', mean(rgb_image,3));
set(gca,'color','k')
caxis([0.0037 0.08])
%caxis([0 1])
colorbar
title('15N12C / (14N12C)')
axis equal, axis tight,% axis off
%

%hold on
%for n = 1:length(posx)
%    text(posx(n)*dx,posy(n)*dx,num2str(n),'color','w')
%end
%plot(posx*dx,posy*dx,'wo');
%for n = 1:length(posx)
%    text(posx(n)*dx,posy(n)*dx,num2str(n),'Color','r')
%end

%% Logarithmic
imagesc(log10(stitchtMat), 'AlphaData', mean(rgb_image,3));
set(gca,'color','k')
%caxis([0.004 0.04])
colormap(jet)
caxis([-3.5 -1])
colorbar
title('15N12C / (14N12C)')

hold on
plot(posx*dx,posy*dx,'wo');
for n = 1:length(posx)
    text(posx(n)*dx,posy(n)*dx,num2str(numb(n)))
end

%%
pNew = p;
pNew.width = size(imNew{1},2);
pNew.height = size(imNew{1},1);
pNew.scale = pNew.width * 0.0781;
p = pNew;
p.im = imNew;
p.planes_aligned = 0;
p.planes = planes;

keep p im m n planes 

%save D:\stitched_grid.mat p -v7.3










