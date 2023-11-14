function    R   =   sure_svt(lambda, sigma, varargin)
    %
    %   SURE_SVT:   Stein's Unbiased Risk Estimate (SURE) for Singular Value
    %               Thresholding (SVT).
    %
    %   USAGE:
    %       R   =   SURE_SVT(LAMBDA, SIGMA, Y)
    %       R   =   SURE_SVT(LAMBDA, SIGMA, S, [M N])    
    %       R   =   SURE_SVT(LAMBDA, SIGMA, S, [M N], IS_REAL)
    %
    %   DESCRIPTION:
    %       R   =   SURE_SVT(LAMBDA, SIGMA, Y)
    %       returns in R the value of SURE for SVT with threshold LAMBDA >
    %       0 evaluated at a matrix Y, which can have real or complex
    %       entries. SIGMA > 0 is the standard deviation for the
    %       noise.
    %       
    %       R   =   SURE_SVT(LAMBDA, SIGMA, S, [M N])
    %       R   =   SURE_SVT(LAMBDA, SIGMA, S, [M N], IS_REAL)
    %       accepts the vector S of observed singular values, and the vector [M N] 
    %       containing the number of rows and columns of the observed matrix. The 
    %       optional flag IS_REAL indicates whether the observed matrix is real-valued 
    %       (if IS_REAL = 1) or not. The default is IS_REAL = 1.
    %
    %   REFERENCES:
    %
    %       Candes, E.J., Sing-Long, C.A., and Trzasko, J.D.,
    %       ``Unbiased Risk Estimates for Singular Value Thresholding''
    %
    %   V1.0:   Oct. 2012.
    %
    if( nargin == 3 ),
        [M N]       =   size(varargin{1});
        is_real     =   isreal(varargin{1});
        [~, s, ~]   =   svd(varargin{1});
        s           =   diag(s(1:min(M,N), 1:min(M,N)));
    elseif( nargin == 4 ),
        s           =   varargin{1};
        M           =   varargin{2}(1);
        N           =   varargin{2}(2);
        is_real     =   1;
    elseif( nargin == 5 ),
        s           =   varargin{1};
        M           =   varargin{2}(1);
        N           =   varargin{2}(2);
        is_real     =   varargin{3};
    end
    
    R           =   sure_div_svt(lambda, s, is_real, [M N]);
    
    %   additional terms
    R           =   -M*N*sigma^2 + sum(min(lambda^2, s.^2)) + 2*sigma^2*R;
    if( ~is_real ),
        R       =   R - M*N*sigma^2;
    end
end
function    x   =   sure_div_svt(lambda, s, is_real, sz)
    svThreshold 	=   1E-8;       %   threshold to determine whether two singular
                                    %   values are the same
    M               =   sz(1);
    N               =   sz(2);
    
    %   *** safeguard
    %       check multiplicities of singular values in a robust manner
    z           =   s(2:end);
    s           =   [s(1) 1];
    Is          =   1;
    while( ~isempty(z) ),
        idx         =   find(abs(z - s(Is, 1)) < svThreshold );
        if( isempty(idx) )
            s           =   [s; [z(1) 1]];
            z(1)        =   [];
            Is          =   Is + 1;
        end
        z(idx)      =   [];
        s(Is, 2)    =   s(Is, 2) + numel(idx);
    end
    clear z
    
    %       warns the user about using SURE with a non-simple, not-full
    %       rank matrix
    if( any( s(:, 1) < svThreshold ) ),
        fprintf('   +   [SURE_SVT] Warning: argument might be rank-deficient.\n')
    end
    if( any( s(:, 2) > 1 ) ),
        fprintf('   +   [SURE_SVT] Warning: argument might have repeated singular values.\n')
    end

    %       find singular values above the threshold
    idx_p   =   ( s(:, 1) > lambda );
    
    if( is_real ),
        x   =   div_svt_real(lambda, s, idx_p, M, N);
    else
        x   =   div_svt_complex(lambda, s, idx_p, M, N);
    end
end
function    x   =   div_svt_real(lambda, s, idx_p, M, N)
    x  	=   0;
    if( any(idx_p)  ),
        x   =   x + sum( 0.5*s(idx_p, 2).*(s(idx_p, 2) + 1) );
        x   =   x + sum( (abs(M-N)*s(idx_p, 2) + 0.5*s(idx_p, 2).*(s(idx_p, 2) - 1)).*(max(0, s(idx_p, 1) - lambda)./s(idx_p, 1)) );
    end
    
    D  	=   zeros(size(s, 1));
    for Ik = 1:size(s, 1),
        D(:, Ik)    =   s(Ik, 2)*s(:, 2).*s(:, 1).*max(0, s(:, 1) - lambda)./(s(:, 1).^2 - s(Ik, 1).^2);
    end
    D( isnan(D) | isinf(D) | abs(D) > 1E6 )     =   0;
    
    x 	=   x + 2*sum(D(:));
end
function    x   =   div_svt_complex(lambda, s, idx_p, M, N)
    x  	=   0;
    if( any(idx_p)  ),
        x   =   x + sum( s(idx_p, 2).^2 );
        x   =   x + sum( (2*abs(M-N) + 1 + s(idx_p, 2).*(s(idx_p, 2) - 1)).*(max(0, s(idx_p, 1) - lambda)./s(idx_p, 1)) );
    end
    
    D  	=   zeros(size(s, 1));
    for Ik = 1:size(s, 1),
        D(:, Ik)    =   s(Ik, 2)*s(:, 2).*s(:, 1).*max(0, s(:, 1) - lambda)./(s(:, 1).^2 - s(Ik, 1).^2);
    end
    D( isnan(D) | isinf(D) | abs(D) > 1E6 )     =   0;
    
    x 	=   x + 4*sum(D(:));
end