classdef AODataCube
    properties
        n
        m
        z
        cube
        thresh
    end
    methods
        function obj = AODataCube(cube)
            obj.cube = cube;
            [m,n,z] = size(cube);
            obj.n = n; % Width
            obj.m = m; % Height
            if m ~= n
                warning('1st and 2nd dimension of datacube are preferably equal')
            end
            obj.z = z;
            obj.thresh = 1e-5;
        end
        function cube2gif(obj,name,delay)
            % This is pretty self explanatory and does most of the work for
            % you. Just call this function from your datacube with a name
            % string and delay in seconds. For example, for a 10fps gif on
            % the datacube 'dc', with a filename of 'animation.gif', you
            % would call:
            %    dc.cube2gif('animation',0.1);
            % And the file 'animation.gif' will appear in whatever
            % directory you run the command from. If you would like to
            % reshape the size of the gif from the defaults, then open
            % figure(1) and adjust to your preference before running the
            % command.
            % For more flexibility (and more work), use the function:
            % frames2gif(...)
            
            figure(1)
            fr = [];
            min_c = min(obj);
            max_c = max(obj);
            for zi = 1:obj.z
                imagesc(obj.cube(:,:,zi))
                set(gca,'YTickLabel','')
                set(gca,'XTickLabel','')
                axis equal
                axis([0.5 obj.n+0.5 0.5 obj.m+0.5])
                caxis([min_c max_c])
                colorbar
                fr = [fr getframe(gcf)];
            end
            frames2gif(fr,name,delay);
        end
        function out = mean(obj)
            out = mean(obj.cube(:));
        end
        function out = max(obj)
            out = max(obj.cube(:));
        end
        function out = min(obj)
            out = min(obj.cube(:));
        end
        function Omega = getOmega(obj)
            phi_hat_a = zeros(obj.n^2*(obj.z-1),1);
            phi_hat_a(:) = obj.cube(:,:,1:obj.z-1);
            tmp = propMatrixBTTB(obj.n,1);
            M = (tmp.m-1):-1:-(tmp.m-1);
            N = (tmp.n-1):-1:-(tmp.n-1);
            Omega = zeros(obj.n^2*(obj.z-1),(2*obj.n-1)^2);
            for mi = 1:(2*obj.n-1)
                for ni = 1:(2*obj.n-1)
                    for zi = 1:(obj.z-1)
                        Omega((obj.n^2*(zi-1)+1):(obj.n^2*zi),(ni-1)*(2*obj.n-1)+mi) = tmp.vsparsemult(M(mi),N(ni),phi_hat_a((obj.n^2*(zi-1)+1):obj.n^2*zi));
                    end
                end
            end
        end
        function Omega = getOffsetOmega(obj,offset)
            phi_hat_a = zeros(obj.n^2*(obj.z-offset),1);
            phi_hat_a(:) = obj.cube(:,:,1:obj.z-offset);
            tmp = propMatrixBTTB(obj.n,1);
            M = (tmp.m-1):-1:-(tmp.m-1);
            N = (tmp.n-1):-1:-(tmp.n-1);
            Omega = zeros(obj.n^2*(obj.z-offset),(2*obj.n-1)^2);
            for mi = 1:(2*obj.n-1)
                for ni = 1:(2*obj.n-1)
                    for zi = 1:(obj.z-offset)
                        Omega((obj.n^2*(zi-1)+1):(obj.n^2*zi),(ni-1)*(2*obj.n-1)+mi) = tmp.vsparsemult(M(mi),N(ni),phi_hat_a((obj.n^2*(zi-1)+1):obj.n^2*zi));
                    end
                end
            end
        end
        function Omega = getTrimmedOffsetOmega(obj,offset)
            Omega = getOffsetOmega(obj,offset);
            Omega(:,obj.getBad) = [];            
        end
        function y = getOffsetY(obj,offset)
            y = zeros(obj.n^2*(obj.z-offset),1);
            y(:) = obj.cube(:,:,(1+offset):obj.z);
        end
        function Omega = getTrimmedOmega(obj)
            Omega = getOmega(obj);
            Omega(:,obj.getBad) = [];
        end
        function bad = getBad(obj)
            bad = (0==getAp(2*obj.n-1));
            bad = bad(:);
        end
        function y = getY(obj)
            y = zeros(obj.n^2*(obj.z-1),1);
            y(:) = obj.cube(:,:,2:obj.z);
        end
    end
end