classdef BCCB < BTTB
    methods
        function obj = BCCB(varargin)
            if nargin==3
                N = varargin{1};
                i = varargin{2};
                j = varargin{3};
                if mod(i,1)~=0 || mod(j,1)~=0
                    error('m and n must be integers');
                end
                mat = zeros(N,N);
                mm = mod(i,N) + 1;
                nn = mod(j,N) + 1;
                mat(mm,nn) = 1;
                mat = [mat(2:end,2:end) mat(2:end,:)
                       mat(:,2:end) mat];
            elseif nargin==1
                defmat = varargin{1};
                [y,x] = size(defmat);
                y = (y+1)/2;
                x = (x+1)/2;
                mat = defmat(1:x,1:x);
                mat(1:x,1:(x-1)) = mat(1:x,1:(x-1)) + defmat(1:x,(x+1):end);
                mat(1:(x-1),1:x) = mat(1:(x-1),1:x) + defmat((x+1):end,1:x);
                mat(1:(x-1),1:(x-1)) = mat(1:(x-1),1:(x-1)) + defmat((x+1):end,(x+1):end);
                mat = padarray(mat,[x-1 x-1],'circular','post');             
            end
            obj@BTTB(mat);
        end
    end
end