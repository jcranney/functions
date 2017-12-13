function fits2gif(fitsfile)
im1 = fitsread(fitsfile);
im = cell(1,length(im1(1,1,:)));
figure(1)
for n = 1:length(im1(1,1,:))
    imagesc(im1(:,:,n))
    axis square
    im{n} = frame2im(getframe(1));
    pause(0.1)
end
filename = [fitsfile(1:end-5) '.gif']; % Specify the output file name
for idx = 1:length(im)
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.5);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.5);
    end
end