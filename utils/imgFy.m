function img = imgFy(y, maskIdxs, dim)
    tmp = zeros(dim(1) * dim(2), dim(3));
    tmp(maskIdxs, :) = y;
    img = reshape(tmp,dim(1),dim(2),[]);
end
