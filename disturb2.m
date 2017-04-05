function d = disturb2(t, varargin)
    n = 2;
    if length(varargin) == 1
        n = varargin{1};
    end
    d(: , 1) = (sin(10*pi*t)+cos(10*pi*t)+sin(9*pi*t)+cos(9*pi*t)+...
        sin(8*pi*t)+cos(8*pi*t)+sin(7*pi*t)+cos(7*pi*t)+...
        sin(6*pi*t)+cos(6*pi*t)+sin(5*pi*t)+cos(5*pi*t)+...
        sin(4*pi*t)+cos(4*pi*t)+sin(3*pi*t)+cos(3*pi*t)-...
        sin(3*pi*t)-cos(3*pi*t)) .* (t>0);
        d(: , 2) = 0.5.*cos(30*t);
        d(: , 3) = 2*sin(10*t);
    d = d(:, 1:n);
end

