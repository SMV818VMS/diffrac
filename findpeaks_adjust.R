##function [p,S,mu]
suppressPackageStartupMessages(library(pracma))

curvefitting=function(x,y,n){ 
# %   POLYFIT Fit polynomial to data. 
# %   POLYFIT(X,Y,N) finds the coefficients of a polynomial P(X) of 
# %   degree N that fits the data, P(X(I))~=Y(I), in a least-squares sense. 
# % 
# %   [P,S] = POLYFIT(X,Y,N) returns the polynomial coefficients P and a 
# %   structure S for use with POLYVAL to obtain error estimates on 
# %   predictions.  If the errors in the data, Y, are independent normal 
# %   with constant variance, POLYVAL will produce error bounds which 
# %   contain at least 50% of the predictions. 
# % 
# %   The structure S contains the Cholesky factor of the Vandermonde 
# %   matrix (R), the degrees of freedom (df), and the norm of the 
# %   residuals (normr) as fields.    
# % 
# %   [P,S,MU] = POLYFIT(X,Y,N) finds the coefficients of a polynomial 
# %   in XHAT = (X-MU(1))/MU(2) where MU(1) = mean(X) and MU(2) = std(X). 
# %   This centering and scaling transformation improves the numerical 
# %   properties of both the polynomial and the fitting algorithm. 
# % 
# %   Warning messages result if N is >= length(X), if X has repeated, or 
# %   nearly repeated, points, or if X might need centering and scaling. 
# % 
# %   See also POLY, POLYVAL, ROOTS. 
# 
# %   Copyright 1984-2002 The MathWorks, Inc. 
# %   $Revision: 5.17 $  $Date: 2002/04/09 00:14:25 $ 
#   
#   % The regression problem is formulated in matrix format as: 
#   % 
# %    y = V*p    or 
# % 
# %          3  2 
# %    y = [x  x  x  1] [p3 
#                        %                      p2 
#                        %                      p1 
#                        %                      p0] 
# % 
# % where the vector p contains the coefficients to be found.  For a 
# % 7th order polynomial, matrix V would be: 
#   % 
# % V = [x.^7 x.^6 x.^5 x.^4 x.^3 x.^2 x ones(size(x))]; 

if(length(x)!=length(y)){ 
  error('X and Y vectors must be the same size.') 
}

##x = x(:) 
##y = y(:) 

mu = c(mean(x), std(x))
x = (x - mu[1])/mu[2]


#% Construct Vandermonde matrix. 
V=matrix(data = 0, nrow = length(x), ncol = n+1)
V[,n+1] = ones(length(x),1)
for (j in n:1){
  V[,j] = x*V[,j+1] 
} 


####% Solve least squares problem, and save the Cholesky factor. 
  QRresult=qr(V,0)

  Q=qr.Q(QRresult)
  R=qr.R(QRresult)

  p=inv(R)%*%(t(Q)%*%y)
  p = c(t(p))##';          % Polynomial coefficients are row vectors by convention. 
 
#        % S is a structure containing three elements: the Cholesky factor of the 
#        % Vandermonde matrix, the degrees of freedom and the norm of the residuals. 
# 
  devuelve=list(p,mu)
  return(devuelve)
}

