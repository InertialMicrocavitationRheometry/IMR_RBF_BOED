
function [idx,D]=knnsearch(varargin)
% KNNSEARCH   Linear k-nearest neighbor (KNN) search
% IDX = knnsearch(Q,R,K) searches the reference data set R (n x d array
% representing n points in a d-dimensional space) to find the k-nearest
% neighbors of each query point represented by eahc row of Q (m x d array).
% The results are stored in the (m x K) index array, IDX. 
%
% IDX = knnsearch(Q,R) takes the default value K=1.
%
% IDX = knnsearch(Q) or IDX = knnsearch(Q,[],K) does the search for R = Q.
%


[Q,R,K,fident] = parseinputs(varargin{:});
% Check outputs
error(nargoutchk(0,2,nargout));
% C2 = sum(C.*C,2)';
[N,M] = size(Q);
L=size(R,1);
idx = zeros(N,K);

D = idx;
if K==1
    % Loop for each query point
    for k=1:N
        d=zeros(L,1);
        for t=1:M
            d=d+(R(:,t)-Q(k,t)).^2;
        end
        if fident
            d(k)=inf;
        end
        [D(k),idx(k)]=min(d);
    end
else
    for k=1:N
        d=zeros(L,1);
        for t=1:M
            d=d+(R(:,t)-Q(k,t)).^2;
        end
        if fident
            d(k)=inf;
        end
        [s,t]=sort(d);
        idx(k,:)=t(1:K);
        D(k,:)=s(1:K);
    end
end
if nargout>1
    D=sqrt(D);
end
function [Q,R,K,fident] = parseinputs(varargin)
% Check input and output
error(nargchk(1,3,nargin));
Q=varargin{1};
if nargin<2
    R=Q;
    fident = true;
else
    fident = false;
    R=varargin{2};
end
if isempty(R)
    fident = true;
    R=Q;
end
if ~fident
    fident = isequal(Q,R);
end
if nargin<3
    K=1;
else
    K=varargin{3};
end
     