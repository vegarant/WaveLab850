function wc = FWT2_CDJV_noP(x,L,N)
% FWT2_CDJV_noP -- 2-d MRA wavelet transform (modified boundary, orthogonal, without preconditioning)
%  Usage
%    wc = FWT2_CDJV_noP(x,L,N)
%  Inputs
%    x     2-d image (n by n array, n dyadic)
%    L     coarse level,  N < 2^L
%    N     degree   integer: 2 or 3 (number of vanishing moments)
%  Outputs
%    wc    2-d wavelet transform
%
%  Description
%    A two-dimensional Wavelet Transform is computed for the
%    array x.  To reconstruct, use IWT2_CDJV.
%
%  See Also
%    IWT2_CDJV
%



	[HPF,LHPEF,RHPEF] = MakeCDJVFilter('HighPass',N);
	[LPF,LLPEF,RLPEF] = MakeCDJVFilter('LowPass',N);
    
    
	[n,J] = quadlength(x);
	wc = x; 
	nc = n;
%     
%     for ix=1:nc,
%         row = wc(ix,1:nc);
%             
%             row(1:N)          =  row(1:N)          * LPREMAT';
%             row(nc:-1:(nc-N+1)) =  row(nc:-1:(nc-N+1)) * RPREMAT';
%             wc(ix,1:nc)=row;
%     end
%     for iy=1:nc,
%         row = wc(1:nc,iy)';
%             
%             row(1:N)          =  row(1:N)          * LPREMAT';
%             row(nc:-1:(nc-N+1)) =  row(nc:-1:(nc-N+1)) * RPREMAT';
%             
%             wc(1:nc,iy)=row';
%     end
%         
    
	for jscal=J-1:-1:L,
		top = (nc/2+1):nc; bot = 1:(nc/2);
		for ix=1:nc,
			row = wc(ix,1:nc);
%             
%             row(1:N)          =  row(1:N)          * LPREMAT';
%             row(nc:-1:(nc-N+1)) =  row(nc:-1:(nc-N+1)) * RPREMAT';
            
			wc(ix,bot) = CDJVDyadDown(row,LPF,LLPEF,RLPEF);
			wc(ix,top) = CDJVDyadDown(row,HPF,LHPEF,RHPEF);
		end
		for iy=1:nc,
			row = wc(1:nc,iy)';
%             
%             row(1:N)          =  row(1:N)          * LPREMAT';
%             row(nc:-1:(nc-N+1)) =  row(nc:-1:(nc-N+1)) * RPREMAT';
%             
            
			wc(top,iy) = CDJVDyadDown(row,HPF,LHPEF,RHPEF)';
			wc(bot,iy) = CDJVDyadDown(row,LPF,LLPEF,RLPEF)'; 
		 end
		nc = nc/2;
    end   
 
% This file was created by Clarice Poon (cmhsp2@cam.ac.uk) 
% This is based on FWT_CDJV.m from Wavelab 850, which was the CDJV forward
% wavelet transform in 1D.

    
