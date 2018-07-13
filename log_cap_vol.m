function d_log_vol = log_cap_vol(n_dim, d_r, d_h)
%% LOG_CAP_VOL Unit hyperspherical cap volume
% Compute the log-volume of a height-d_h, radius-d_r cap on the unit
% n_dim-ball.
% 
% Input: n_dim = dimension, natural number
%        d_r = radius of ball
%        d_h = cap height, in [0,d_r]
%
% Output: d_log_vol = log(volume of cap)

%% Strategy
% Vol = 1/2 * ...
%       pi^(n_dim/2) * ...
%       1/gamma(1+n_dim/2) * ...
%       d_r^n_dim * ...
%       I_{(2*d_r*d_h-d_h^2)/d_r^2}((n_dim+1)/2, 1/2)
% where I_x(a,b) is the regularized incomplete beta function. 
% So:
% log(Vol) = log(1/2) + ...
%            (n_dim/2)*log(pi) - ...
%            gammaln(1+n_dim/2) + ...
%            n_dim*log(d_r) + ...
%            log(betainc(d_x,d_a,d_b)).
% where:
% d_x = (2*d_r*d_h-d_h^2)/d_r^2
% d_a = (n_dim+1)/2
% d_b = 1/2
% betainc(.) = MATLAB Regularized incomplete beta function.

%% Substituting log(betainc(.))
% log(betainc(.)) may not be precise enough. It may be necessary to replace
% this term with more precise ones:
% log(betainc(.)) = log(B_(d_x)(d_a,d_b)/beta(d_a,d_b))
%                 = log(B_(d_x)(d_a,d_b)) - betaln(d_a,d_b)
%                 = log(2F1(d_a+d_b, 1; d_a+1; d_x)) + ...
%                   d_a*log(d_x) + ...
%                   d_b*log(1-d_x) - ...
%                   log(d_a) - ...
%                   betaln(d_a,d_b).
% 2F1( . , . ; . ; . ) = Gauss' hypergeometric function.
    
    d_log_vol_sphere = (n_dim/2)*log(pi) - ...
                       gammaln(1+n_dim/2) + ...
                       n_dim*log(d_r);
    d_log_vol = [];
    if(d_h < d_r)
        d_x = (2*d_r*d_h-d_h^2)/d_r^2;
        d_a = (n_dim+1)/2;
        d_b = 1/2;
        log_b_x = log(betainc(d_x,d_a,d_b));
        
        d_log_vol = log(1/2) + ...
            d_log_vol_sphere + ...
            log_b_x;
    elseif(d_h == d_r)
        d_log_vol = d_log_vol_sphere-log(2);
    else
        d_log_vol = NaN;
    end
    return;
end
