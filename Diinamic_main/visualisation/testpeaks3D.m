function tests = testpeaks3D (image,p,fit)
% function tests = testpeaks3D (image,p,fit)
% ChiTst:(1-alpha)-quantile of the fit
% Ftest
% Marianne Renner avril 09 for SPTrack v4.0
%----------------------------------------------------------

n = prod(size(image));

% chi^2- test of the fit
  residues = (image - fit);
  %residues = residues ./ sqrt(noise^2+fit);
  residues = residues ./ sqrt(fit);
  ChiTst   = (n-1) * std(residues(:))^2;
  ChiTst   = 1 - chipf (n-1,ChiTst);
 
 % F-test for the residues 
  FTest = std(residues(:))^2;
  if FTest<1, FTest=1/FTest; end
  FTest = (n-1) / (2*(n-1)*FTest);
  FTest = 2 * betainc (FTest,0.5*(n-1),0.5*(n-1));
  if FTest>1, FTest=2-FTest; end

  tests=[ChiTst,FTest];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ChiProbFun = chipf(n,z)
nt = 1000;
t = 0:z/nt:z;
ChiProbFun=1/2^(n/2)/gamma(n/2)*sum(t.^(n/2-1).*exp(-t/2))*z/nt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

