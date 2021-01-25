function orientim = orient_img_calc(normim ,gradientsigma, blocksigma, orientsmoothsigma)

%gradientsigma=10; blocksigma=15; orientsmoothsigma=15;
%gradientsigma=1; blocksigma=5; orientsmoothsigma=5;
[rows,cols] = size(normim);
% Calculating image gradients
sze = fix(6*gradientsigma);
if ~mod(sze,2);
    sze = sze+1;
end
f = fspecial('gaussian', sze, gradientsigma); % Generate Gaussian filter.
[fx,fy] = gradient(f);                        % Gradient of Gausian.
Gx = filter2(fx, normim); % Gradient of the image in x
Gy = filter2(fy, normim); % and y
% Estimate the local orientation of each block
%

%
Gxx = Gx.^2;
Gxy = Gx.*Gy;
Gyy = Gy.^2;
sze = fix(6*blocksigma);
if ~mod(sze,2);
    sze = sze+1;
end
f = fspecial('gaussian', sze, blocksigma);
%f = ones(40,40);

Gxx = filter2(f, Gxx);  % Smoothed
Gxy = 2*filter2(f, Gxy);
Gyy = filter2(f, Gyy);
denom = sqrt(Gxy.^2 + (Gxx - Gyy).^2);
sin2theta = Gxy./denom;            % Sine and cosine
cos2theta = (Gxx-Gyy)./denom;
sze = fix(6*orientsmoothsigma);
if ~mod(sze,2);
    sze = sze+1;
end
f = fspecial('gaussian', sze, orientsmoothsigma);
%f = ones(40,40);

cos2theta = filter2(f, cos2theta); % Smoothed sine and cosine
sin2theta = filter2(f, sin2theta);
orientim = mod(atan2(-sin2theta,cos2theta),2*pi)/2;
