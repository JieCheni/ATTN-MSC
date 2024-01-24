function TU = mode_n1n2_unfold(varargin)
% mode_n1n2_unfold unfold a tensor into a tensor of size In1xIn2x...IR.
%     Note that IR is the product of remaining dimensions,IR = ΠIi,i∈[1,T.ndims] and
%     i≠n1,...,nn
%     TU = mode_n1n2_unfold(T) unfold tensor T into a matrix by arranging the first 
%     mode as row and other modes as column in little-endian order, just
%     like mode_n_unfolding().
%
%     TU = mode_n1n2_unfold(T,n) unfold tensor T into a tensor unfolded of size
%     In1xIn2x...xInnxIR in little-endian order,n is a vector. Note that IR is the
%     product of remaining dimensions,IR = ΠIi,i∈[1,T.ndims] and
%     i≠n1,...,nn
%
%     TU = mode_n1n2_unfold(T,n1,n2...nn) unfold tensor T into a tensor unfolded of size
%     In1xIn2x...xInnxIR in little-endian order. Note that IR is the
%     product of remaining dimensions,IR = ΠIi,i∈[1,T.ndims] and
%     i≠n1,...,nn
% Examples
%         T = tensor(rand(2,4,6,8));
%         MT = mode_n1n2_unfold(T);
%         MT = mode_n1n2_unfold(T,2);
%         MT = mode_n1n2_unfold(T,[3,2]);
%         MT = mode_n1n2_unfold(T,3,2);

    % Input must be Tensor class
    if nargin == 0
        error("Input tensor must not be empty!")
    end
    if isa(varargin{1},'tensor')
        TU = varargin{1}.data;
    elseif isnumeric(varargin{1})
        TU = varargin{1};
    else
        error("Input tensor must be tensor class data or high dimension matrix!")
    end

    sz = size(TU);
    ndim = length(sz);

    if nargin == 1
        modeN = 1;
    elseif nargin == 2 
        if all(isnumeric(varargin{2})) && all(varargin{2} == int8(varargin{2})) && all(min(varargin{2})>0) && all(max(varargin{2})<=ndim)
            modeN = int8(varargin{2});
        else
            error("Please input correct mode vector(integer between 1 and ndim)!")
        end
    else
        % nargin >= 3
        if isnumeric([varargin{2:end}]) && all([varargin{2:end}] == int8([varargin{2:end}])) && all(min([varargin{2:end}])>0) && all(max([varargin{2:end}])<=ndim)
            modeN = int8([varargin{2:end}]);
        else
            error("Please input correct mode vector(integer between 1 and ndim)!")
        end
    end

    % Get the last mode IR
    idx = true(ndim,1);
    idx(modeN) = false;
    mode_R = 1:ndim;
    mode_R = mode_R(idx);
    
    % Litte endian order
    TU = permute(TU,[modeN,mode_R]);
    mode_sz_N = sz(modeN);
    mode_sz_R = prod(sz(mode_R));
    modeList = [mode_sz_N,mode_sz_R];
    TU = reshape(TU,modeList);
   
end