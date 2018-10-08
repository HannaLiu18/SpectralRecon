function ProjData = EllipseProj(Phantom, structProj)

N_Ellipse = size(Phantom,1);
ProjData = zeros(structProj.nAngle,structProj.nChannel);

fAngle = 2 * pi/structProj.nAngle;
fDetWidth = 2 / structProj.nChannel;

for ii = 1:N_Ellipse
    A = Phantom(ii,1);           % Amplitude change for this ellipse
    asq = Phantom(ii,2)^2;       % a^2
    bsq = Phantom(ii,3)^2;       % b^2
    x0 = Phantom(ii,4);          % x offset
    y0 = Phantom(ii,5);          % y offset
    phi = Phantom(ii,6)*pi/180;  % rotation angle in radians
    
    for jj = 1:structProj.nAngle
        theta = (jj-1) * fAngle;
        thetaphi = theta - phi;
        cospt = cos(thetaphi); sinpt = sin(thetaphi);
        TruncA =  x0 * cos(phi) + y0 * sin(phi);
        TruncB = -x0 * sin(phi) + y0 * cos(phi);
        
        for kk = 1:structProj.nChannel
            b = fDetWidth * (kk - structProj.nChnCenter);
            
            TempA = (b * cospt - TruncA);
            TempB = (b * sinpt - TruncB);
            AA = sinpt * sinpt / asq + cospt * cospt / bsq;
            BB = 2 * sinpt * TempA / asq - 2 * cospt * TempB / bsq;
            CC = TempA * TempA / asq + TempB * TempB / bsq - 1;
            
            Delta = BB^2 - 4 * AA * CC;
            if (Delta > 0)
                ProjData(jj,kk) = ProjData(jj,kk) + A * sqrt(Delta) /(2*AA);
            end
        end
    end
end

ProjData = ProjData * structProj.fFOV;
end