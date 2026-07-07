// pchip.rs (splines library)

//! Pchip interpolation method (TODO)

use num::{Float};

use crate::spline::{Spline};
use crate::{diff};

pub fn pchip<T: Float + std::fmt::Debug>(xx: &[T], yy: &[T])->Spline<T>{
	let hh : Vec<T> = diff(xx).collect();
	let delta : Vec<T> = diff(yy).zip(hh.iter()).map(|(dy,dx)| dy/ *dx).collect();
	let ss = slopes_pchip(xx, yy, &delta);
	return Spline::new(xx, yy, &ss);
}

/*
% Check and adjust input data
[x,y,sizey] = chckxy(x,y);

% Compute slopes
h = diff(x);
m = prod(sizey);
delta = diff(y,1,2)./repmat(h,m,1);
slopes = zeros(size(y));
for r = 1:m
    if isreal(delta)
        slopes(r,:) = pchipSlopes(x,y(r,:),delta(r,:));
    else
        realSlopes = pchipSlopes(x,y(r,:),real(delta(r,:)));
        imagSlopes = pchipSlopes(x,y(r,:),imag(delta(r,:)));
        slopes(r,:) = complex(realSlopes,imagSlopes);
    end
end
*/

/// TODO
pub fn slopes_pchip<T: Float>(xx: &[T], yy: &[T], delta: &[T])->Vec<T>{
	todo!();
}

// Matlab code

/*
function d = pchipSlopes(x,y,del)
% Derivative values for PCHIP.
% d = pchipslopes(x,y,del) computes the first derivatives, d(k) = P'(x(k)).

%  Special case n = 2, use linear interpolation.
n = length(x);
if n == 2
    d = repmat(del(1),size(y));
    return
end

%  Slopes at interior points.
%  d(k) = weighted average of del(k-1) and del(k) when they have the same sign.
%  d(k) = 0 when del(k-1) and del(k) have opposites signs or either is zero.

k = find(sign(del(1:n-2)).*sign(del(2:n-1)) > 0);

h = diff(x);
hs = h(k)+h(k+1);
w1 = (h(k)+hs)./(3*hs);
w2 = (hs+h(k+1))./(3*hs);
dmax = max(abs(del(k)), abs(del(k+1)));
dmin = min(abs(del(k)), abs(del(k+1)));

d = zeros(size(y));
d(k+1) = dmin./conj(w1.*(del(k)./dmax) + w2.*(del(k+1)./dmax));

%  Slopes at end points.
%  Set d(1) and d(n) via non-centered, shape-preserving three-point formulae.

d(1) = ((2*h(1)+h(2))*del(1) - h(1)*del(2))/(h(1)+h(2));
if sign(d(1)) ~= sign(del(1))
    d(1) = 0;
elseif (sign(del(1)) ~= sign(del(2))) && (abs(d(1)) > abs(3*del(1)))
    d(1) = 3*del(1);
end

d(n) = ((2*h(n-1)+h(n-2))*del(n-1) - h(n-1)*del(n-2))/(h(n-1)+h(n-2));
if sign(d(n)) ~= sign(del(n-1))
    d(n) = 0;
elseif (sign(del(n-1)) ~= sign(del(n-2))) && (abs(d(n)) > abs(3*del(n-1)))
    d(n) = 3*del(n-1);
end

*/