function [ A ] = polyDef( x,y,z,pOrd )

% Check input arguments
if nargin == 4
    dim = 3;
elseif nargin == 3
    pOrd = z;   
    dim = 2;
else    
    error('Not enough input arguments');   
end
    
    
switch dim
    
    % 2 dimensional polynomial
    case 2
        
        % Order 0
        if pOrd == 0
            A = zeros(size(x,1),1);
            A(:,1) = ones(size(x));
        end
        
        % Order 1
        if pOrd == 1
            A = zeros(size(x,1),3);
            A(:,1) = ones(size(x));
            A(:,2) = (x.^1) .* (y.^0);
            A(:,3) = (x.^0) .* (y.^1);
            
        % Order 2
        elseif pOrd == 2
            A = zeros(size(x,1),6);
            A(:,1)  = ones(size(x));
            A(:,2)  = (x.^1).*(y.^0);
            A(:,3)  = (x.^0).*(y.^1);
            A(:,4)  = (x.^2).*(y.^0);
            A(:,5)  = (x.^1).*(y.^1);
            A(:,6)  = (x.^0).*(y.^2);
            
        % Order 3    
        elseif pOrd == 3
            A = zeros(size(x,1),10);
            A(:,1)  = ones(size(x));
            A(:,2)  = (x.^1).*(y.^0);
            A(:,3)  = (x.^0).*(y.^1);
            A(:,4)  = (x.^2).*(y.^0);
            A(:,5)  = (x.^1).*(y.^1);
            A(:,6)  = (x.^0).*(y.^2);
            A(:,7)  = (x.^3).*(y.^0);
            A(:,8)  = (x.^2).*(y.^1);
            A(:,9)  = (x.^1).*(y.^2);
            A(:,10) = (x.^0).*(y.^3);
            
        % Order 4    
        elseif pOrd == 4
            A = zeros(size(x,1),15);
            A(:,1)  = ones(size(x));
            A(:,2)  = (x.^1).*(y.^0);
            A(:,3)  = (x.^0).*(y.^1);
            A(:,4)  = (x.^2).*(y.^0);
            A(:,5)  = (x.^1).*(y.^1);
            A(:,6)  = (x.^0).*(y.^2);
            A(:,7)  = (x.^3).*(y.^0);
            A(:,8)  = (x.^2).*(y.^1);
            A(:,9)  = (x.^1).*(y.^2);
            A(:,10) = (x.^0).*(y.^3);
            A(:,11) = (x.^4).*(y.^0);
            A(:,12) = (x.^3).*(y.^1);
            A(:,13) = (x.^2).*(y.^2);
            A(:,14) = (x.^1).*(y.^3);
            A(:,15) = (x.^0).*(y.^4);
            
        % Order 5    
        elseif pOrd == 5
            A = zeros(size(x,1),21);
            A(:,1)  = ones(size(x));
            A(:,2)  = (x.^1).*(y.^0);
            A(:,3)  = (x.^0).*(y.^1);
            A(:,4)  = (x.^2).*(y.^0);
            A(:,5)  = (x.^1).*(y.^1);
            A(:,6)  = (x.^0).*(y.^2);
            A(:,7)  = (x.^3).*(y.^0);
            A(:,8)  = (x.^2).*(y.^1);
            A(:,9)  = (x.^1).*(y.^2);
            A(:,10) = (x.^0).*(y.^3);
            A(:,11) = (x.^4).*(y.^0);
            A(:,12) = (x.^3).*(y.^1);
            A(:,13) = (x.^2).*(y.^2);
            A(:,14) = (x.^1).*(y.^3);
            A(:,15) = (x.^0).*(y.^4);
            A(:,16) = (x.^5).*(y.^0);
            A(:,17) = (x.^4).*(y.^1);
            A(:,18) = (x.^3).*(y.^2);
            A(:,19) = (x.^2).*(y.^3);
            A(:,20) = (x.^1).*(y.^4);
            A(:,21) = (x.^0).*(y.^5);
        else
            error('Only a polynomial of up to 5th order is supported for 2D polynomial');
        end
     
    % 3 dimensional polynomal    
    case 3
        
       % Order 0
        if pOrd == 0
            A = zeros(size(x,1),1);
            A(:,1)  = ones(size(x));
        
        % Order 1
        elseif pOrd == 1
            A = zeros(size(x,1),3);
            A(:,1) = ones(size(x));
            A(:,2) = (x.^1) .* (y.^0) .* (z.^0);
            A(:,3) = (x.^0) .* (y.^1) .* (z.^0);
            A(:,4) = (x.^0) .* (y.^0) .* (z.^1);
            
        % Order 2
        elseif pOrd == 2
            A = zeros(size(x,1),10);
            A(:,1) = ones(size(x));
            A(:,2) = x;
            A(:,3) = y;
            A(:,4) = z;
            A(:,5) = x.^2;
            A(:,6) = x.*y;
            A(:,7) = y.^2;
            A(:,8) = x.*z;
            A(:,9) = y.*z;
            A(:,10) = z.^2;
            
        % Order 3    
        elseif pOrd == 3
            A = zeros(size(x,1),20);
            A(:,1) = ones(size(x));
            A(:,2) = x;
            A(:,3) = y;
            A(:,4) = z;
            A(:,5) = x.^2;
            A(:,6) = x.*y;
            A(:,7) = y.^2;
            A(:,8) = x.*z;
            A(:,9) = y.*z;
            A(:,10) = z.^2;
            A(:,11) = x.^3;
            A(:,12) = (x.^2).*y;
            A(:,13) = x.*(y.^2);
            A(:,14) = y.^3;
            A(:,15) = (x.^2).*z;
            A(:,16) = x.*y.*z;
            A(:,17) = (y.^2).*z;
            A(:,18) = x.*(z.^2);
            A(:,19) = y.*(z.^2);
            A(:,20) = z.^3;
            
        % Order 4    
        elseif pOrd == 4
            A = zeros(size(x,1),35);
            A(:,1) = ones(size(x));
            A(:,2) = x;
            A(:,3) = y;
            A(:,4) = z;
            A(:,5) = x.^2;
            A(:,6) = x.*y;
            A(:,7) = y.^2;
            A(:,8) = x.*z;
            A(:,9) = y.*z;
            A(:,10) = z.^2;
            A(:,11) = x.^3;
            A(:,12) = (x.^2).*y;
            A(:,13) = x.*(y.^2);
            A(:,14) = y.^3;
            A(:,15) = (x.^2).*z;
            A(:,16) = x.*y.*z;
            A(:,17) = (y.^2).*z;
            A(:,18) = x.*(z.^2);
            A(:,19) = y.*(z.^2);
            A(:,20) = z.^3;
            A(:,21) = x.^4;
            A(:,22) = (x.^3).*y;
            A(:,23) = (x.^2).*(y.^2);
            A(:,24) = x.*(y.^3);
            A(:,25) = y.^4;
            A(:,26) = (x.^3).*z;
            A(:,27) = (x.^2).*y.*z;
            A(:,28) = x.*(y.^2).*z;
            A(:,29) = (y.^3).*z;
            A(:,30) = (x.^2).*(z.^2);
            A(:,31) = x.*y.*(z.^2);
            A(:,32) = (y.^2).*(z.^2);
            A(:,33) = x.*(z.^3);
            A(:,34) = y.*(z.^3);
            A(:,35) = z.^4;

        else
            error('Only a polynomial of up to 4th order is supported for 3D polynomial');
        end
end
end
        

