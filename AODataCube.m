classdef AODataCube
    properties
        n
        z
        cube
        thresh
    end
    methods
        function obj = AODataCube(cube)
            obj.cube = cube;
            obj.n = size(cube,1);
            if size(cube,1) ~= size(cube,2)
                error('1st and 2nd dimension of datacube must be equal')
            end
            obj.z = size(cube,3);
            obj.thresh = 1e-5;
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