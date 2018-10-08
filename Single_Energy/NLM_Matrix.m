function MaskM = NLM_Matrix(Nim,Image,sigma,NPatch,NNeighbour)
    XLoc = [];
    YLoc = [];
    Data = [];
    nn = 0;Nx = Nim; Ny = Nim;   
    for xx = (1 + NPatch):(Nx - NPatch)
        % xx
        for yy = (1 + NPatch):(Nx - NPatch)
            % yy
            Patch_Ref = Image(xx-NPatch:xx + NPatch,yy-NPatch:yy + NPatch);
            nLoc_Ref = (yy-1) * Nx + xx;
            WeightSum = 0;
            NNStart = nn;
            for xxPatch = max(xx - NNeighbour,1 + NPatch):min((xx + NNeighbour),Nx - NPatch)
                for yyPatch = max(yy - NNeighbour,1 + NPatch) :min((yy + NNeighbour),Ny - NPatch)
                    nLoc_Comp = (yyPatch-1) * Nx + xxPatch;
                    
                    if (nLoc_Comp ~= nLoc_Ref)
                    Patch_Comp = Image(xxPatch-NPatch:xxPatch + NPatch,...
                        yyPatch-NPatch:yyPatch + NPatch);
                    Weight = exp(-norm(Patch_Ref(:) - Patch_Comp(:))^2 / sigma);
%                     if (Weight < 10^(-10))
%                         disp('Something unexpected!');
%                         pause;
%                     en=[
                    Weight = sqrt(Weight);
                    if (Weight < 0.01)
                        Weight = 0.01;
                    end
                    if (Weight > 0.5)
                        Weight = 0.5;
                    end
                    WeightSum = WeightSum + Weight;
                    XLoc(nn * 2 + 1:nn * 2 + 2) = [nn + 1,nn + 1];
                    YLoc(nn * 2 + 1:nn * 2 + 2) = [nLoc_Ref, nLoc_Comp];
                    Data(nn * 2 + 1:nn * 2 + 2) = [Weight,-Weight];
                    nn = nn + 1;
                    end
                end
            end 
            Data(NNStart * 2 + 1:nn * 2) = Data(NNStart * 2 + 1:nn * 2) / WeightSum;
%             % To set a possible smallest value for weights
%             Data(NNStart * 2 + 1:2:nn * 2) = max( 0.02, Data(NNStart * 2 + 1:2:nn * 2));
%             Data(NNStart * 2 + 2:2:nn * 2) = min(-0.02, Data(NNStart * 2 + 1:2:nn * 2));
%             WeightSum = sum(Data(NNStart * 2 + 1:2:nn * 2));
%             Data(NNStart * 2 + 1:nn * 2) = Data(NNStart * 2 + 1:nn * 2) / WeightSum;
        end
    end
    nn
    MaskM = (sparse(XLoc,YLoc,Data,nn, Nx * Ny, nn * 2));
end