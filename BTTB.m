classdef BTTB
    properties
        DefMat
        m
        n
    end
    methods
        function obj = BTTB(varargin)
            if nargin == 1
                mat = varargin{1};
                obj.DefMat = mat;
                obj.m = (size(mat,2)+1)/2;
                obj.n = (size(mat,1)+1)/2;
                if mod(obj.m,1)~= 0 || mod(obj.n,1)~=0
                    error('Incompatible defining matrix dimensions for BTTB')
                end
            elseif nargin == 2
                obj.m = varargin{1};
                obj.n = varargin{2};
                obj.DefMat = zeros(2*obj.m-1,2*obj.n-1);
            end
        end
        function A = getFullMatrix(obj)
            Acell = cell(obj.n,obj.n);
            nidx = (obj.n-1):-1:-(obj.n-1);
            blks = toeplitz(0:-1:-(obj.n-1),0:(obj.n-1));
            for ni = 1:(2*obj.n-1)
                a_ = obj.DefMat(ni,:);
                a = fliplr(a_(1:obj.m));
                b = a_(obj.m:end);
                tmp = toeplitz(b,a);
                [Acell{(blks==nidx(ni))}] = deal(tmp);               
            end
            A = cell2mat(Acell);
        end
        function [sz,varargout] = size(obj,varargin)
            sz_ = [obj.n*obj.m obj.n*obj.m];
            if nargin == 2
                idx = varargin{1};
                sz = sz_(idx);
            else
                sz = sz_;
            end
            if nargout == 2
                sz = sz_(1);
                varargout = {sz_(2)};
            end
        end         
        function z = getElement(obj,idx)
            if length(idx) ~= 2
                [tmpa, tmpb] = sub2ind(size(obj),idx);
                idx = [tmpa tmpb];
            end
            iq = floor((idx(1)-1)/obj.m);
            jq = floor((idx(2)-1)/obj.m);
            ioff = mod((idx(1)-1),obj.m);
            joff = mod((idx(2)-1),obj.m);
            z = obj.DefMat(obj.n+iq-jq,obj.m+ioff-joff);
        end
        function C = mtimes(A,B)
            %% BTTB * 
            if isa(A,'BTTB') && isa(B,'double')
                if all(size(B)==1)
                    C = BTTB(A.DefMat*B);
                    return
                elseif size(A,2)~=size(B,1)
                    error('Internal dimensions must agree');
                end
                C = zeros(size(A,1),size(B,2));
                for k = 1:A.m
                    for l = 1:A.n
                        a = A.getACol(k,l);
                        b = A.getBRow(B,k,l);
                        C = C + kron(a,b);                
                    end
                end
            elseif isa(B,'BTTB') && isa(A,'double')
                if all(size(A)==1)
                    C = BTTB(B.DefMat*A);
                    return
                elseif size(A,2)~=size(B,1)
                    error('Internal dimensions must agree');
                end
                C = sparse(size(B,2),size(A,1));
                tmp = A';
                A = B';
                B = tmp;
                for k = 1:A.m
                    for l = 1:A.n
                        a = A.getACol(k,l);
                        b = A.getBRow(B,k,l);
                        C = C + kron(a,b);                
                    end
                end
                C = C';
            elseif isa(A,'BTTB') && isa(B,'BTTB')
                if size(A,2)~=size(B,1)
                    error('Internal dimensions must agree');
                end
                B = B';
                C = zeros(size(A));
                for k = 1:A.m
                    for l = 1:A.n
                        a = A.getACol(k,l);
                        b = B.getACol(k,l)';
                        C = C + kron(a,b);                
                    end
                end
            end
        end
        function C = ctranspose(obj)
            C = BTTB(rot90(obj.DefMat,2));
        end
        function col = getACol(obj,k,l)
            col = zeros(obj.m*obj.n,1);
            for mi = 1:obj.m
                for ni = 1:obj.n
                    col((mi-1)*obj.m+ni) = col((mi-1)*obj.m+ni) + obj.DefMat(-(k-mi)+obj.m,-(l-ni)+obj.n);
                end
            end
        end
        function b = getBRow(obj,B,k,l)
            b = B((k-1)*obj.m+l,:);
        end
        function C = plus(A,B)
            if isa(A,'BTTB')
                C = A.getFullMatrix + B;
            elseif isa(B,'BTTB')
                C = A + B.getFullMatrix;
            else
                error('Neither A nor B are BTTB');
            end
        end
        function C = minus(A,B)
            C = A+(-1.0*B);
        end
        function obj = assignLambdas(obj,m,n,lambda)
            if length(m)~=length(n) || length(m)~=length(lambda)
                error('When using vector inputs, all must be same length')
            end
            for mi = 1:length(m)
                obj.DefMat(obj.m + m(mi),obj.n + n(mi)) = lambda(mi);
            end
        end
        function C = mpower(A,B)
            if ~(isa(A,'BTTB') && isa(B,'double') && length(B)==1) || mod(B,1)~=0
                error('Incompatible exponent');
            else
                C = speye(size(A));
            end
            for ni = 1:B
                C = C*A;
            end
        end
        function A = vsparsemult(obj,alpha,beta,phi)
            A = zeros(obj.m*obj.n,1);
            alpha = -alpha;
            beta = -beta;
            mi = (max(alpha,0)+1):(max(alpha,0)+obj.m-abs(alpha));
            ni = (max(beta,0)+1):(max(beta,0)+obj.n-abs(beta));
            [Mi,Ni] = ndgrid(mi,ni);
            idxold = sub2ind([obj.m obj.n],Mi(:),Ni(:)); % these previous two lines could be done faster
            mi = (max(-alpha,0)+1):(max(-alpha,0)+obj.m-abs(alpha));
            ni = (max(-beta,0)+1):(max(-beta,0)+obj.n-abs(beta));
            [Mi,Ni] = ndgrid(mi,ni);
            idxnew = sub2ind([obj.m obj.n],Mi(:),Ni(:));
            A(idxnew) = phi(idxold);
        end
        function v = eig(obj)
            v = eig(obj.getFullMatrix);
        end
    end
end
