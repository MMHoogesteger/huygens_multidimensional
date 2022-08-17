function [ak_n] = normalizeFourierCoefficients(ak,an)
R = size(ak,2);
ejtau = -1i*abs(an)/an;
ak_n = zeros(size(ak));
for r = 1:R
    ak_n(:,r) = ejtau^r*ak(:,r);
end

end

