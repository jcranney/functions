function frames2gif(fr,filename,varargin)
% FRAMES2GIF Takes frames in fr and writes them to [filename '.gif'] in the
% current working directory. The third argument is the delay between frames
% in the final gif (default 0.5s).
%   FRAMES2GIF(FR,FILENAME) writes FR to FILENAME with delay of 0.5
%   FRAMES2GIF(FR,FILENAME,DELAY) writes FR to FILENAME with DELAY
%
% Jesse Cranney, 2016
if nargin == 3
    delay = varargin{1};
else
    delay = 0.5;
end
im = cell(length(fr),1);
for n = 1:length(fr)
    im{n} = frame2im(fr(n));
end
filename = [filename '.gif']; % Specify the output file name
for idx = 1:length(im)
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',delay);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',delay);
    end
end