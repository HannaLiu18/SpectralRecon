function MaskM = FirstDerivativeM(Nim)
XLoc = [];
YLoc = [];
Data = [];
nn = 0;Nx = Nim; Ny = Nim;
for yy = 1:(Ny-1)
    for xx = 1:(Nx-1)
        nLoc = (yy-1) * Nx + xx;
        XLoc(nn * 4 + 1:nn * 4 + 4) = [nn *2 + 1,nn*2 + 1,nn*2+2,nn*2+2];
        YLoc(nn * 4 + 1:nn * 4 + 4) = [nLoc, nLoc + 1, nLoc, nLoc + Nx];
        Data(nn * 4 + 1:nn * 4 + 4) = [1,-1,1,-1];
        nn = nn + 1;
    end
end
XLine = nn * 2;
nn = nn * 4;

yy = Ny;
for xx = 1:(Nx-1)
    nLoc = (yy-1) * Nx + xx;
    XLoc(nn + 1:nn + 2) = [XLine + 1, XLine + 1];
    YLoc(nn + 1:nn + 2) = [nLoc, nLoc + 1];
    Data(nn + 1:nn + 2) = [1,-1];
    nn = nn + 2;
    XLine = XLine + 1;
end

xx = Nx;
for yy = 1:(Ny-1)
    nLoc = (yy-1) * Nx + xx;
    XLoc(nn + 1:nn + 2) = [XLine + 1, XLine + 1];
    YLoc(nn + 1:nn + 2) = [nLoc, nLoc + Nx];
    Data(nn + 1:nn + 2) = [1,-1];
    nn = nn + 2;
    XLine = XLine + 1;
end

MaskM = (sparse(XLoc,YLoc,Data,XLine, Nx * Ny, nn));
end