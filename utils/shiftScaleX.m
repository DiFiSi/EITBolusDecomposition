function Xout = shiftScaleX(X,mu,scale)
    Xout = (X - mu) ./ scale;
end
