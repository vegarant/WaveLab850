function x= IWT_CDJV_noP(wc,L,N)
% IWT_CDJV_noP -- Inverse Wavelet Transform  (boundary corrected, without postconditioning)
%  Usage
%    x = IWT_CDJV_noP(wc,L,N)
%  Inputs
%    wc   1-d wavelet transform
%    L    Level of V_0;  L << J
%    N    Degree of Daubechies Filters
%  Outputs
%    x    1-d signal: length(y) = 2^J
%

	[HPF,LHPEF,RHPEF] = MakeCDJVFilter('HighPass',N);
	[LPF,LLPEF,RLPEF] = MakeCDJVFilter('LowPass',N);
%
    wcoef = ShapeAsRow(wc);
	[~,J] = dyadlength(wcoef) ;
	beta = wcoef(1:(2^(L))); 
	for j=L:(J-1),
	   alfa = CDJVDyadUp(wcoef(dyad(j)),HPF,LHPEF,RHPEF); 
	   beta = CDJVDyadUp(beta,LPF,LLPEF,RLPEF) + alfa;
	end
	x = beta;
% 	x(1:N) = (beta(1:N)) *  LPOSTMAT';
% 	x(n:-1:(n-N+1)) = beta(n:-1:(n-N+1)) * RPOSTMAT';
%
    x = ShapeLike(x,wc);
    
% This file is IWT_CDJV.m from Wavelab 850, but without preconditioning.
