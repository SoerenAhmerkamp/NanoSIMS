%% Simple detection of FISH positive cells

load alignedExternal.mat
cd mat
figure
subplot(1,2,1)
imshow(im./255)
axis equal, axis tight
%enrich_12C15N = sum(double(p.im{4})./double(p.im{3}),3);

load 12C15N.mat
v12C15N = IM;
load 12C14N.mat
v12C14N = IM;

enrich_12C15N = v12C15N ./ v12C14N;

subplot(1,2,2)
imagesc(enrich_12C15N)
axis equal, axis tight
%%
dapi = im;
dapi(:,:,2) = 0;
dapi(:,:,1) = 0;

cells = im(:,:,3).*im(:,:,1)./im(:,:,2)./255;
cells(cells>5) = 0;
imagesc(cells)
caxis([0.2 1])
axis equal, axis tight

mask = cells>0.5;
mask = double(mask);
%mask(mask == 0) = 0.01;
imagesc(mask);


[B,L] = bwboundaries(mask,'noholes');
%imshow(label2rgb(L, @parula, [.5 .5 .5]))
celleri15N_spec = [];
hold on
for k = 1:length(B)
   boundary = B{k};
   %plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 0.5)
   temp = L;
   temp(~(L==k)) = 0;
   temp(temp>0) = 1;
   cell_sp = enrich_12C15N.*temp;
   celleri15N_spec(end+1) = mean(cell_sp(cell_sp>0));
end
%%
plant = im(:,:,2)./im(:,:,3)./255;
plant(plant>0.02)=0;
imagesc(plant)
maskplant = plant>0.01;
maskplant = double(maskplant);
imagesc(maskplant)

%%
subplot(1,3,1)
mask(mask==0) = NaN;
celleri15N = medfilt2(enrich_12C15N.*mask,[4 4]);
imagesc(celleri15N)
caxis([0 100])

subplot(1,3,2)
maskplant(maskplant==0) = NaN;
plant15N = medfilt2(enrich_12C15N.*maskplant,[3 3]);
imagesc(plant15N)
caxis([0 100])

subplot(1,3,3)
imshow(im./255)%
%[maskBG] = drawrectangle();
maskBG = ~(im(:,:,1) & im(:,:,2) & im(:,:,3));
BG15N = medfilt2(enrich_12C15N.*maskBG,[10 10]);
imagesc(BG15N)
caxis([0 100])

figure
celi = padarray(celleri15N_spec',length(plant15N(:))-length(celleri15N_spec),NaN,'post');

boxdata = [celi plant15N(:) BG15N(:)];
boxplot(boxdata)
ylim([0.0 0.12])