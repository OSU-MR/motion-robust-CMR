function [ v_1, v_0 ] = post_v(rho_b, rho_v, tauSq_b, tauSq_v, sigmaSq, gamma )
%POST_V Calculates the posterior probabilities for the hidden support
%variable v
%   Detailed explanation goes here

sigmaSq_p = tauSq_b + tauSq_v + sigmaSq;

% % Loop through the index
% for ind = 1:length(rho_b)
%
%     % Calculate the unnormalized probability of no velocity
%     v_0(ind) = (1- gamma(ind))./(pi*sigmaSq_p(ind)).*exp(-abs(rho_v(ind)-rho_b(ind)).^2./sigmaSq_p(ind));
%
%     % Use different equation depending on the size of the input of the
%     % bessel function
%     if (2*abs(rho_b(ind)).*abs(rho_v(ind)))/sigmaSq_p(ind) < 10
%
%         % Exact equation for small values
%         v_1(ind) = gamma(ind)/(pi*sigmaSq_p(ind))*exp(-(abs(rho_b(ind))^2 + abs(rho_v(ind))^2)/sigmaSq_p(ind));
%         b(ind) = besseli(0,2*abs(rho_b(ind))*abs(rho_v(ind))/sigmaSq_p(ind));
%         v_1(ind) = v_1(ind).*b(ind);
%     else
%         % Approximation for larg Values
%         v_1(ind) = gamma(ind)/(2*sqrt(pi^3.*sigmaSq_p(ind)*abs(rho_b(ind))*abs(rho_v(ind))))...
%             *exp(-(abs(rho_b(ind))^2 + abs(rho_v(ind))^2 - 2*abs(rho_b(ind))*abs(rho_v(ind)) )/sigmaSq_p(ind));
%     end
% end

v_0 = (1- gamma)./(pi*sigmaSq_p).*exp(-abs(rho_v-rho_b).^2./sigmaSq_p);

v_1 = gamma./(2*sqrt(pi^3.*sigmaSq_p.*abs(rho_b).*abs(rho_v)))...
    .*exp(-(abs(rho_b).^2 + abs(rho_v).^2 - 2*abs(rho_b).*abs(rho_v) )./sigmaSq_p);

% v_1 = gamma./(pi*sigmaSq_p).*exp(-(abs(rho_b).^2 + abs(rho_v).^2)./sigmaSq_p);
% b = besseli(0,2*abs(rho_b).*abs(rho_v)./sigmaSq_p);
% v_1 = v_1.*b;

% Normalize
sums = v_0 + v_1 +eps;
v_0 = v_0./sums + eps;
v_1 = v_1./sums + eps;

end

