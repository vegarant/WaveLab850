function x = IWT2_CDJV_noP(wc,L,N)
% IWT2_CDJV_noP -- Inverse 2-d MRA wavelet transform (modified boundary, orthogonal, without postconditioning)
%  Usage
%    x = IWT2_CDJV_noP(wc,L,N)
%  Inputs
%    wc    2-d wavelet transform [n by n array, n dyadic]
%    L     coarse level
%    N     degree   integer: 2 or 3 (number of vanishing moments)
%  Outputs
%    x     2-d signal reconstructed from wc
%
%
%
%

    [HPF,LHPEF,RHPEF] = MakeCDJVFilter('HighPass',N);
	[LPF,LLPEF,RLPEF] = MakeCDJVFilter('LowPass',N);
%
    

	[~,J] = quadlength(wc);
	x = wc; 
	nc = 2^(L+1);
    for jscal=L:J-1,
        top = (nc/2+1):nc; bot = 1:(nc/2); all = 1:nc;
		for iy=1:nc,
            x_interm =   CDJVDyadUp(x(bot,iy)',LPF,LLPEF,RLPEF)  ...
					   + CDJVDyadUp(x(top,iy)',HPF,LHPEF,RHPEF); 
            x(all,iy) = x_interm';
            
		end
		for ix=1:nc,
            x_interm =   CDJVDyadUp(x(ix,bot),LPF,LLPEF,RLPEF)  ...
					   + CDJVDyadUp(x(ix,top),HPF,LHPEF,RHPEF); 
            x(ix,all) = x_interm;
            
		end
		nc = 2*nc;
    end
	
% This file was created by Clarice Poon (cmhsp2@cam.ac.uk) 
% This is based on IWT_CDJV.m from Wavelab 850, which was the CDJV inverse
% wavelet transform in 1D with preconditioning.


