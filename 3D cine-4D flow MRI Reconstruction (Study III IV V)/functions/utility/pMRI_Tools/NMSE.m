function [nmse] = NMSE(xtrue,xhat)

nmse = 20*log10(norm(xhat(:)-xtrue(:))/norm(xtrue(:)));
end