classdef geomGrad < handle
    properties (SetAccess = public)
        validLenslet;
        resolution;
        pupil;
        edgePix;
        
    end
    properties (SetObservable=true)
        % measurements
        slopes=0;
        % gradient matrix
        gradOp;
    end
    methods
        %% Constructor
        function obj = geomGrad(validLenslet, pupil, varargin)
            %error(nargchk(1,3,nargin));
            p = inputParser;
            p.addRequired('validLenslet', @(x) isnumeric(x) || islogical(x));
            p.addRequired('pupil', @(x) isnumeric(x) || islogical(x));
            p.addOptional('edgePix', 0, @isnumeric);
            p.parse(validLenslet, pupil,varargin{:});
            
            obj.validLenslet = validLenslet;
            obj.resolution  = size(pupil,1);
            obj.pupil = pupil;
            obj.edgePix = p.Results.edgePix;
            
            computeGradOp(obj);
        end
        
        %% Relay
        % Method that allows compatibility with the overloaded mtimes
        % operator, allowing things like source=(source.*tel)*wfs
        function relay(obj, src)
            obj.slopes = obj.gradOp*src.phase(obj.pupil);
        end
        
        %% Compute gradients
        function obj = computeGradOp(obj)
            
            nSlopes   = 2*nnz(obj.validLenslet);
            G         = zeros(nSlopes,nnz(obj.pupil)); iL = 1;
            nLenslet  = size(obj.validLenslet,1);
            nPx_subap = obj.resolution/nLenslet;
            nPx_subap2= nPx_subap^2;
            for j = 1:nLenslet
                for i = 1:nLenslet
                    % for a valid lenslet
                    if obj.validLenslet(i,j) == 1
                        
                        % valid sub-aperture pixels
                        subapMask = zeros(obj.resolution);
                        subapMask((i-1)*nPx_subap+1:i*nPx_subap,(j-1)*nPx_subap+1:j*nPx_subap) = 1;
                        subapMask = subapMask.*obj.pupil;
                        
                        % X gradients
                        
                        maskLet = zeros(obj.resolution); nPixUsed = 0;
                        if nnz(subapMask) == nPx_subap2 % fully illuminated sub-aperture
                            maskLet((i-1)*nPx_subap+1:i*nPx_subap,(j-1)*nPx_subap+1) = -1;
                            maskLet((i-1)*nPx_subap+1:i*nPx_subap,j*nPx_subap) = 1;
                            if obj.edgePix
                                nPixUsed = nPx_subap2;
                            else
                                nPixUsed = nPx_subap2-nPx_subap;
                            end
                        else % partially illuminated sub-aperture
                            for ik = (i-1)*nPx_subap+1:i*nPx_subap
                                lig = subapMask(ik,:);
                                if nnz(lig) >= 2
                                    jl = find(lig,1,'first');
                                    jr = find(lig,1,'last');
                                    nPixUsed = nPixUsed + (jr-jl+obj.edgePix);
                                    maskLet(ik,jl) = -1;% left pixel
                                    maskLet(ik,jr) = 1;% right pixel
                                end
                            end
                            
                        end
                        maskLet = maskLet*nPx_subap/nPixUsed/pi;
                        % rasterize
                        M = maskLet(:);
                        % keep only valid pixels
                        Mpup = M(obj.pupil,1);
                        % populate gradient matrix
                        if obj.edgePix
                            G(iL,:) = Mpup'*(nPx_subap/(nPx_subap-1));
                        else
                            G(iL,:) = Mpup';
                        end
                        
                        % Y gradients
                        maskLet = zeros(obj.resolution); nPixUsed = 0;
                        if nnz(subapMask) == nPx_subap2 % fully illuminated sub-aperture
                            maskLet((i-1)*nPx_subap+1,(j-1)*nPx_subap+1:j*nPx_subap) = -1;
                            maskLet(i*nPx_subap, (j-1)*nPx_subap+1:j*nPx_subap) = 1;
                            if obj.edgePix
                                nPixUsed = nPx_subap2;
                            else
                                nPixUsed = nPx_subap2-nPx_subap;
                            end
                        else
                            for jk = (j-1)*nPx_subap+1:j*nPx_subap
                                col = subapMask(:,jk);
                                if nnz(col) >= 2
                                    ih = find(col,1,'first');
                                    il = find(col,1,'last');
                                    nPixUsed = nPixUsed + (il-ih+obj.edgePix);
                                    maskLet(ih,jk) = -1;% pixel du HAUT
                                    maskLet(il,jk) = 1;% pixel du BAS
                                end
                            end
                        end
                        maskLet = maskLet*nPx_subap/nPixUsed/pi;
                        % rasterize
                        M = maskLet(:);
                        % keep only valid pixels
                        Mpup = M(obj.pupil,1);
                        % populate gradient matrix
                        if obj.edgePix
                            G(iL+nSlopes/2,:) = Mpup'*(nPx_subap/(nPx_subap-1));
                        else
                            G(iL+nSlopes/2,:) = Mpup';
                        end
                        iL = iL+1;
                        
                    end
                end
            end
            [i,j,s] = find(G);
            obj.gradOp = sparse(i,j,s,nSlopes,nnz(obj.pupil));
        end
    end
    
end
