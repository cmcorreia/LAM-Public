classdef pyramid < handle
    % PYRAMID Creates a pyramid WFS object
    %
    % pwfs = pyramid(nLenslet,nPx) creates a pyramid WFS object equivalent
    % to a (nLenslet X nLenslet) Shack-Hartmann WFS with an incoming
    % wavefront with a (nPx X nPx) resolution
    %
    % pwfs = pyramid(...,'modulation',m) modulates the pyramid with an
    % amplitude m lambda/D
    %
    % pwfs = pyramid(...,'binning',n) bins the detector by factor n before
    % computing the slopes
    
    % pwfs = pyramid(...,'obstructionRatio',r) creates a
    % validIntensityPupil with a central obstruction with ratio r
    
    properties
        
        nLenslet                    % # of lenslet
        resolution;                 % resolution of the single quadrant detector images
        camera;                     % the detector used at the end of the chain
        referenceSlopesMap;         % the slopes map of reference
        framePixelThreshold = -inf; % frame pixel threshold
        obstructionRatio;           % central obstrction ratio [0...1]
        slopesDisplayHandle;        % slopes display handle
        slopesListener;             % slope listener
        intensityDisplayHandle      % intensity display handle
        tag = 'PYRAMID';            % tag
    end
    properties(GetAccess = 'public', SetAccess = 'private')
        validActuator;              % map of validActuators [default at corners of sub-apertures]
        lightMap;                   % map of the light passed through the pyramid
        slopesMap;                  % the slopes as read by processing the sensor
        wave;                       % incoming wave, a square complex matrix
        validSlopes;                % logical insindex indicating the validSlopes
        validIntensityPupil;        % logical index indicating the validIntensityPupil
        validLenslet;               % logical index coresponding to valid pixel points but 
                                    % over dimensions of nLensletxnLenslet
                                    % (for compatability with shackHartmann
                                    % functions)
        pyrMask;                    % pyramid face transmitance and phase
        slopesUnits=1;              % normalisation factor to have the 1-to-1 gain btw input and output
    end
    
    properties(SetAccess=private,SetObservable=true)
        slopes;                     % the slopes as read by processing the sensor
    end
    
    properties (Dependent)
        modulation;                 % amplitude of the modulation in \lambda/D units, a decimal number
        alpha;                      % angle of incoming rays
        c;                          % Multiplicative term of the Nyquist sampling (which is 2 by default). With  c = 2 (default) the quadrants are spread by a factor 2. End result :: c=2 -> 4pixels across \lambda/D
        binning;                    % binning factor, default 1
    end
    
    properties (Dependent,SetAccess=private)
        nSlope;                     % # of slopes
    end
    
    properties (Access=private)

        p_modulation;
        p_alpha;
        p_c;
        p_binning = 1;              % binning factor, default 1
        isInitialized;
        u;
        pxSide;
        o;
        r;
        nTheta;
        q;
        I4Q4;                       % intensity map on the 4 pyramid quadrants before detector binning
        log;
    end
    
    %%
    methods
        %% Constructor
        function pwfs = pyramid(nLenslet,resolution,varargin)
            p = inputParser;
            addRequired(p,'nLenslet',@isnumeric);
            addRequired(p,'resolution',@isnumeric);
            addParameter(p,'modulation',0,@isnumeric);
            addParameter(p,'binning',1,@isnumeric);
            addParameter(p,'c',2,@isnumeric);
            addParameter(p,'alpha',pi/2,@isnumeric);
            addParameter(p,'obstructionRatio',0,@isnumeric);
            parse(p,nLenslet,resolution,varargin{:})
            
            pwfs.nLenslet           = p.Results.nLenslet;
            pwfs.resolution         = p.Results.resolution;
            pwfs.p_alpha            = p.Results.alpha;
            pwfs.obstructionRatio   = p.Results.obstructionRatio;
            pwfs.c                  = p.Results.c;
            pwfs.modulation         = p.Results.modulation;
            pwfs.camera             = detector(2*pwfs.c*nLenslet);
            pwfs.binning            = p.Results.binning;
            pwfs.validSlopes        = [pwfs.validIntensityPupil pwfs.validIntensityPupil];
            
            pwfs.camera.frameGrabber=pwfs; %
            
            pwfs.referenceSlopesMap = zeros(pwfs.nLenslet*pwfs.c/pwfs.binning,pwfs.nLenslet*2*pwfs.c/pwfs.binning);
            
            pwfs.slopesListener     = addlistener(pwfs,'slopes','PostSet',...
                @(src,evnt) pwfs.slopesDisplay );
            pwfs.slopesListener.Enabled = false;
            pwfs.log                = logBook.checkIn(pwfs);
        end
        
        %% Destructor
        function delete(pwfs)
            if isvalid(pwfs.camera)
                delete(pwfs.camera)
            end
            checkout(pwfs.log,pwfs);
        end
        
        %% SETS AND GETS
        %% Get nSlope
        function out = get.nSlope(pwfs)
            out = sum(pwfs.validSlopes(:));
        end
        
        %% Get and Set modulation
        function out = get.modulation(pwfs)
            out = pwfs.p_modulation;
        end
        function set.modulation(pwfs,val)
            pwfs.p_modulation = val;
            [uu,vv] = ndgrid((0:(pwfs.pxSide-1))./pwfs.pxSide);
            [pwfs.o,pwfs.r] = cart2pol(uu,vv);
            pwfs.nTheta = round(2*pi*pwfs.c*pwfs.modulation);
            setModulationData(pwfs,1)
        end
        
        function setModulationData(pwfs,nWave)
            pwfs.q   = zeros(pwfs.pxSide,pwfs.pxSide,nWave);
            pwfs.I4Q4 = zeros(pwfs.pxSide,pwfs.pxSide,...
                nWave,pwfs.nTheta); %
        end
        %% Get and Set validIntensityPupil
        function setValidIntensityPupil(pwfs,val)
            if nargin == 2
                pwfs.validIntensityPupil = logical(val);
            else
                pwfs.validIntensityPupil = utilities.piston((pwfs.nLenslet)/pwfs.p_binning,...
                    pwfs.nLenslet*pwfs.c/pwfs.p_binning,'type','logical','xOffset',0,'yOffset',0);
                if pwfs.obstructionRatio
                    centralHole = utilities.piston(floor(pwfs.obstructionRatio*(pwfs.nLenslet)/pwfs.p_binning),...
                        pwfs.nLenslet*pwfs.c/pwfs.p_binning,'type','logical','xOffset',0,'yOffset',0);
                else
                    centralHole = 0;
                end
                pwfs.validIntensityPupil = logical(pwfs.validIntensityPupil - centralHole);
            end
        end
        
        %% Get and Set validLenslet
        function setValidLenslet(pwfs)
           
            % Sets indices validLenslet for compatibility with other
            % functions for shackHartmann etc.  This gives the valid
            % lenslet/pixel pupil over a grid of nLenslet X nLenslet
            % dimensions.
            pwfs.validLenslet = utilities.piston(pwfs.nLenslet/pwfs.p_binning,...
                'type','logical','xOffset',0,'yOffset',0);
            if pwfs.obstructionRatio
                centralHole = utilities.piston(floor(pwfs.obstructionRatio*pwfs.nLenslet/pwfs.p_binning),...
                    pwfs.nLenslet/pwfs.p_binning,'type','logical','xOffset',0,'yOffset',0)
            else
                centralHole = 0;
            end
            pwfs.validLenslet = logical(pwfs.validLenslet - centralHole);
            
        end
        %% Get and Set alpha
        function out = get.alpha(pwfs)
            out = pwfs.p_alpha;
        end
        function set.alpha(pwfs,val)
            pwfs.p_alpha = val;
            makePyrMask(pwfs);
        end
        
        %% Get and Set c
        function out = get.c(pwfs)
            out = pwfs.p_c;
        end
        function set.c(pwfs,val)
            pwfs.p_c = val;
            pwfs.u = 1+pwfs.resolution*(2*pwfs.c-1)/2:pwfs.resolution*(2*pwfs.c+1)/2;
            pwfs.pxSide = pwfs.resolution*2*pwfs.p_c;
            makePyrMask(pwfs);
            setValidIntensityPupil(pwfs);
            setValidLenslet(pwfs);
            pwfs.validSlopes = [pwfs.validIntensityPupil pwfs.validIntensityPupil];
        end
        
        %% Get and Set binning
        function out = get.binning(pwfs)
            out = pwfs.p_binning;
        end
        function set.binning(pwfs,val)
            pwfs.p_binning = val;
            pwfs.isInitialized = false;
            setValidIntensityPupil(pwfs);
            pwfs.validSlopes = [pwfs.validIntensityPupil pwfs.validIntensityPupil];
         end

        %% Get valid actuators
        function val = get.validActuator(obj)
            nElements            = 2*obj.nLenslet+1; % Linear number of lenslet+actuator
            validLensletActuator = zeros(nElements);
            index                = 2:2:nElements; % Lenslet index
            idx0 = floor(obj.nLenslet*(obj.c-1)/2);
            w = idx0+1:idx0+obj.nLenslet;
            validLensletActuator(index,index) = obj.validIntensityPupil(w,w);
            for xLenslet = index
                for yLenslet = index
                    if validLensletActuator(xLenslet,yLenslet)==1
                        xActuatorIndice = [xLenslet-1,xLenslet-1,...
                            xLenslet+1,xLenslet+1];
                        yActuatorIndice = [yLenslet-1,yLenslet+1,...
                            yLenslet+1,yLenslet-1];
                        validLensletActuator(xActuatorIndice,yActuatorIndice) = 1;
                    end
                end
            end
            index = 1:2:nElements; % Actuator index
            val   = logical(validLensletActuator(index,index));
        end
       
        %% Relay
        % Method that allows compatibility with the overloaded mtimes
        % operator, allowing things like source=(source.*tel)*wfs
        function relay(pwfs, src)
            %% RELAY pyramid to source relay
            pyramidTransform(pwfs,src.catWave);
            grab(pwfs.camera,pwfs.lightMap);
            dataProcessing(pwfs);
        end
        
        function varargout = slopesDisplay(pwfs,varargin)
            %% SLOPESDISPLAY pyramide slopes display
            %
            % slopesDisplay(obj) displays image plot of the slopes
            %
            % slopesDisplay(obj,'PropertyName',PropertyValue) displays
            % image plot of the slopes and set the properties of the
            % graphics object quiver
            %
            % h = slopesDisplay(obj,...) returns the graphics handle
            n = pwfs.nLenslet*pwfs.c;
            data = reshape( pwfs.slopesMap , n , 2*n, []);
            data = bsxfun( @times, data, pwfs.validSlopes);
            data = reshape( data, n , []);
            if ishandle(pwfs.slopesDisplayHandle)
                set(pwfs.slopesDisplayHandle,...
                    'Cdata',data,varargin{:})
            else
                pwfs.slopesDisplayHandle = imagesc(data,varargin{:});
                ax = gca;
                pos = get(ax,'position');
                axis xy equal tight
                ylabel(colorbar('location','EastOutside'),'Pixel')
                set(ax,'position',pos)
            end
            if nargout>1
                varargout{1} = pwfs.slopesDisplayHandle;
            end
        end
        
        function varargout = intensityDisplay(pwfs,varargin)
            %% INTENSITYDISPLAY pyramid lenslet intensity display
            %
            % intensityDisplay(obj) displays the intensity of the lenslets
            %
            % intensityDisplay(obj,'PropertyName',PropertyValue) displays
            % the intensity of the lenslets and set the properties of the
            % graphics object imagesc
            %
            % h = intensityDisplay(obj,...) returns the graphics handle
            %
            % See also: imagesc
            n = size(pwfs.camera.frame,1);
            frame = pwfs.camera.frame(1:n/2,:) + pwfs.camera.frame(1+n/2:n,:);
            frame = reshape(frame',size(frame,1)*2,[]);
            n = size(frame,1);
            intensity = frame(1:n/2,:) + frame(1+n/2:n,:);
            if ishandle(pwfs.intensityDisplayHandle)
                set(pwfs.intensityDisplayHandle,...
                    'CData',intensity,...
                    varargin{:})
            else
                pwfs.intensityDisplayHandle = ...
                    imagesc(intensity,varargin{:});
                ax = gca;
                pos = get(ax,'position');
                axis equal tight xy
                colorbar('location','EastOutside')
                set(ax,'position',pos)
            end
            if nargout>0
                varargout{1} = obj.intensityDisplayHandle;
            end
        end
        
        %% INITIALISATION :: referecence slopes and gain calibration
        function INIT(pwfs)
            %% INIT pyramid inialization
            %
            % pwfs.INIT sets the reference slopes, maps and gains
            
            pwfs.referenceSlopesMap = pwfs.slopesMap + pwfs.referenceSlopesMap;
            pwfs.slopesUnits = 1;
            pwfs.slopesMap = pwfs.slopesMap*0;
            pwfs.slopes = pwfs.slopes*0;
            gainCalibration(pwfs)
            pwfs.isInitialized = true;
        end
        
        %This function updates the detector with whatever light there
        %is in the pyramid, and then procedes to process the output of
        %the detector in order to obtain the slopes
        function uplus(pwfs)
            %% UPLUS + Update operator
            %
            % +obj grabs a frame and computes the slopes
            %
            % obj = +obj returns the pyramid object
            
            grab(pwfs.camera,pwfs.lightMap);
            dataProcessing(pwfs);
            
        end
        
        %This function computes the slopes corresponding to a flat
        %incoming light wave.
        function dataProcessing(pwfs)
            %% DATAPROCESSING Processing a SH-WFS detector frame
            %
            % dataProcessing(obj) computes the WFS slopes
            
            n = size(pwfs.camera.frame,1)/pwfs.binning;
            frame = pwfs.camera.frame;
            if isfinite(pwfs.framePixelThreshold)
                frame = frame - pwfs.framePixelThreshold;
            end
            frame(frame<0) = 0;
            I4Q = utilities.binning(frame,size(pwfs.camera.frame)/pwfs.binning);             % binning
            I4Q = reshape(I4Q,n,n,[]);
            
            n = size(I4Q,1);
            nFrames = size(I4Q,3);
            if nFrames > 1
                im = mat2cell(I4Q,[n/2 n/2],[n/2 n/2],nFrames);
            else
                im = mat2cell(I4Q,[n/2 n/2],[n/2 n/2]);
            end
            computeSlopes(pwfs, im{1}, im{2}, im{4}, im{3});% modif
        end
        
        %% gain calibration
        function gainCalibration(pwfs)
            nPx = pwfs.resolution;
            ngs = source('wavelength',photometry.R);
            tel = telescope(1,'resolution',nPx);
            zer = zernike(3,'resolution', nPx);
            sx = zeros(1,5);
            sy = zeros(1,5);
            for i = 1:5
                zer.c = (i-3)*0.1/(ngs.waveNumber);
                ngs = ngs.*tel*zer*pwfs;
                sx(i) = mean(pwfs.slopes(1:end/2));
                sy(i) = mean(pwfs.slopes(end/2+1:end));
            end
            Ox_in = 4*(-2:2)*0.1;%/ngs.waveNumber;
            Ox_out = sy;
            slopesLinCoef = polyfit(Ox_in,Ox_out,1);
            pwfs.slopesUnits = 1/slopesLinCoef(1);
            
            ngs = ngs.*tel*pwfs;
        end
        %% Hilbert transform approximation for the slopes
        function hilbertSlopes(pwfs,telescope,source)
            P=telescope.pupil;
            phi=angle(source.wave);
            %            Sx = ???PHx (P??) + P??Hx (P) ??? Hxy (P)Hy (P??) + Hxy (P??)Hy (P),
            %            Sy = ???PHy (P??) + P??Hy (P) ??? Hxy (P)Hx (P??) + Hxy (P??)Hx (P).
            %            Sx = -P*imag(hilbert(P*phi)) + P*phi*imag(hilbert(P)) - imag(hilbert(imag(hilbert(transpose(P)))))*imag(hilbert(transpose(P*phi))) + imag(hilbert(imag(hilbert(transpose(P*phi)))))*imag(hilbert(transpose(P))) ;
            %            Sy = -P*imag(hilbert(transpose(P*phi))) + P*phi*imag(hilbert(transpose(P))) - imag(hilbert(imag(hilbert(transpose(P)))))*imag(hilbert(P*phi)) + imag(hilbert(imag(hilbert(transpose(P*phi)))))*imag(hilbert(P)) ;
            %Sx = -P.*imag(hilbert(P.*phi)) + P.*phi.*imag(hilbert(P)) - imag(hilbert(transpose(imag(hilbert(P))))).*imag(hilbert(transpose(P.*phi))) + imag(hilbert(transpose(imag(hilbert(P.*phi))))).*imag(hilbert(transpose(P))) ;
            %Sy = -P*imag(hilbert(transpose(P*phi))) + P*phi*imag(hilbert(transpose(P))) - imag(hilbert(imag(hilbert(transpose(P)))))*imag(hilbert(P*phi)) + imag(hilbert(imag(hilbert(transpose(P*phi)))))*imag(hilbert(P)) ;
            myHilbertX = @(x) imag(hilbert(x));
            myHilbertY = @(x) imag(hilbert(x')');
            myHilbertXY = @(x) imag(hilbert(imag(hilbert(x')')));
            Sx = -P.*myHilbertX(P.*phi) + P.*phi.*myHilbertX(P) - myHilbertXY(P).*myHilbertY(P.*phi) + myHilbertXY(P.*phi).*myHilbertY(P);
            Sy = -P.*myHilbertY(P.*phi) + P.*phi.*myHilbertY(P) - myHilbertXY(P).*myHilbertX(P.*phi) + myHilbertXY(P.*phi).*myHilbertX(P) ;
            pwfs.slopes=[Sx,Sy];
            
        end
        
    end
    
    methods (Access=private)
        
        %% Propagate the wavefront transformed by the pyram WFS to the detector
        function pwfs = pyramidTransform(pwfs,wave)
            [n1,n2,n3] = size(wave);
            nWave = n2*n3/n1;
            if size(pwfs.q,3)~=nWave
                setModulationData(pwfs,nWave)
            end
            
            pwfs.q(pwfs.u,pwfs.u,:) = reshape(wave, n1,[],nWave);
            
            if pwfs.modulation>0
                for kTheta = 1:pwfs.nTheta
                    theta = (kTheta-1)*2*pi/pwfs.nTheta;
                    fftPhasor = exp(-1i.*pi.*4*pwfs.modulation*...
                        pwfs.c*pwfs.r.*cos(pwfs.o+theta));
                    buf = bsxfun(@times,pwfs.q,fftPhasor);
                    buf = bsxfun(@times,fft2(buf),pwfs.pyrMask);
                    pwfs.I4Q4(:,:,:,kTheta) = abs(fft2(buf)).^2;
                    
                end
                %                 figure
                %                 imagesc(tmp)
                %                 axis square
                %                 title(sprintf('Modulation=%d;c=%d',pwfs.modulation,pwfs.c))
                I4Q = sum(pwfs.I4Q4,4);
                
            else
                I4Q = bsxfun(@times,fft2(pwfs.q),pwfs.pyrMask);
                I4Q = abs(fft2(I4Q)).^2;
            end
            pwfs.lightMap = reshape(I4Q,pwfs.pxSide,...
                pwfs.pxSide*nWave);
        end
        
        %% make the pyramid phase mask
        function makePyrMask(pwfs)
            
            [fx,fy] = freqspace(pwfs.pxSide,'meshgrid');
            fx = fx.*floor(pwfs.pxSide/2);
            fy = fy.*floor(pwfs.pxSide/2);
            %pym = zeros(pxSide);
            
            % pyramid face transmitance and phase for fx>=0 & fy>=0
            mask  = heaviside(fx).*heaviside(fy);
            phase = -pwfs.alpha.*(fx+fy);
            pym   = mask.*exp(1i.*phase);
            % pyramid face transmitance and phase for fx>=0 & fy<=0
            mask  = heaviside(fx).*heaviside(-fy);
            phase = -pwfs.alpha.*(fx-fy);
            pym   = pym + mask.*exp(1i.*phase);
            % pyramid face transmitance and phase for fx<=0 & fy<=0
            mask  = heaviside(-fx).*heaviside(-fy);
            phase = pwfs.alpha.*(fx+fy);
            pym   = pym + mask.*exp(1i.*phase);
            % pyramid face transmitance and phase for fx<=0 & fy>=0
            mask  = heaviside(-fx).*heaviside(fy);
            phase = -pwfs.alpha.*(-fx+fy);
            pym   = pym + mask.*exp(1i.*phase);
            pwfs.pyrMask   = fftshift(pym./sum(abs(pym(:))));
        end
        
        %% Compute slopes from 4 intensity maps
        function computeSlopes(pwfs, I1, I2, I3, I4)
            
            % normalisation options
            %    1) normalisation pixel-wise by the intensity
            I = (I1+I2+I3+I4);      %
            %    2) normalisation by the integrated flux
            I2D = utilities.toggleFrame(I,2);
            I = sum(I2D(pwfs.validIntensityPupil))*ones(size(I));
            
            SyMap = (I1-I2+I4-I3)./I;
            SxMap = (I1-I4+I2-I3)./I;
            
            pwfs.slopesMap = bsxfun(@minus,[SxMap,SyMap],pwfs.referenceSlopesMap);
            
            [n1,n2,n3] = size(pwfs.slopesMap);
            slopesMap_ = reshape(pwfs.slopesMap,n1*n2,n3);
            pwfs.slopes = slopesMap_(pwfs.validSlopes(:),:)*pwfs.slopesUnits;
            pwfs.slopesMap = reshape( slopesMap_ , n1, n2*n3);
        end
        
    end
end




















