function laggedX = lagmatrix(x,lags)

laggedX=[];
for i=1:lags+1   
   laggedX = [laggedX x(i:end-lags+i-1,:)];
end