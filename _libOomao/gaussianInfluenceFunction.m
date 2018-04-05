classdef gaussianInfluenceFunction < handle
    % INFLUENCEFUNCTION Create an influence function object
    %
    % obj = influenceFunction('monotonic',mech) creates a cubic Bezier
    % influence function monotically decreasing from 1 to 0 with the
    % mechanical coupling mech
    %
    % obj = influenceFunction('overshoot',mech) creates a cubic Bezier
    % influence function with a negative overshoot and the mechanical
    % coupling mech
    %
    % Try show(obj) to view a cut of the influence function
    
    properties
        % mechanicalCoupling
        mechCoupling;
        % spline polynomials coefs
        splineP;
        % path listener
        gaussianListener;
        % modes
        modes
        % actuators coordinates
        actuatorCoord;
        % condition for influence function centring
        influenceCentre;
        % influence function tag
        tag = 'GAUSSIAN INFLUENCE FUN';
    end
    
    properties (SetObservable=true)
        % path
        gaussian;
        pitch; % actuator pitch in physical units
    end
    
    
    properties (Access=private)
        % influence function display handle;
        displayHandle
        log
    end
    
    methods
        
        %% Constructor
        function obj = gaussianInfluenceFunction(mech,pitch, influenceCentre)
            if ~exist('pitch','var')
                pitch = 1; % default pitch value
            end
            if ~exist('influenceCentre','var')
                influenceCentre = 0;
            end
            
            obj.mechCoupling = mech;
            obj.pitch = pitch;
            obj.influenceCentre = influenceCentre;
            obj.gaussianListener = addlistener(obj,'gaussian','PostSet',...
                @(src,evnt) obj.show );
            obj.gaussianListener.Enabled = false;
            obj.log = logBook.checkIn(obj);
            obj.splineP = mech;
            
        end
        
        %% Destructor
        function delete(obj)
            if ishandle(obj.displayHandle)
                delete(get(obj.displayHandle,'parent'));
            end
            checkOut(obj.log,obj)
        end
        
        function display(obj)
            %% DISP Display object information
            %
            % disp(obj) prints information about the influenceFunction
            % object
            fprintf('___ %s ___\n',obj.tag)
            fprintf('  . mechanical coupling: %3.1f\n',...
                obj.mechCoupling);  
            fprintf('  . pitch in meters: %3.1f\n',...
                obj.pitch);  
            fprintf('----------------------------------------------------\n')
        end
                
        %% Set spline
        function set.splineP(obj,val)
            
            % Notes: only using simple setup: x axis is linearly defined,
            % unlike in bezier (influenceFunction) class.  May change this
            % later.
            
            c = 1/sqrt(log(1/val));
            df = 1e-10;
            mx = sqrt(-log(df)*c^2);
            x = linspace(-mx,mx,1001);
            f = exp(-x.^2/c^2);
            
            obj.gaussian(:,1) = x*obj.pitch;
            obj.gaussian(:,2) = f;
            obj.splineP = spline(x*obj.pitch,f);
        end
        
        function out = mtimes(obj,c)
            %% MTIMES Influence function multiplication
            %
            % v = obj*c multiply the influence functions in obj by the
            % vector c
            
            out = obj.modes*c;
        end
        function out = mldivide(obj,c)
            %% MLDIVIDE Influence function projection
            %
            % v = obj\c project the influence functions onto the vector c
            
            out = obj.modes\c;
        end
        
        function show(obj,varargin)
            %% SHOW 
            
            %figure
                plot(obj.gaussian(:,1),obj.gaussian(:,2))
                hold on
                plot(obj.pitch*ones(100,1), linspace(0,1,100),'k--')
                plot(-obj.pitch*ones(100,1), linspace(0,1,100),'k--')
                plot(obj.gaussian(:,1), obj.mechCoupling*ones(size(obj.gaussian(:,1))),'k--')
                axis tight
                xlabel('position [m]')
                ylabel('normalised influence function')
            
        end
        
        function setInfluenceFunction(obj,nIF,resolution,validActuator,ratioTelDm,offset, diameter)
            %% SETINFLUENCEFUNCTION
            
            %             if nargin<5
            %                 ratioTelDm = 1;
            %             end
            %             z = linspace(-1,1,nIF)*(nIF-1)/2;
            if isempty(obj.actuatorCoord)
                
                if obj.influenceCentre==0
                    xIF = linspace(-1,1,nIF)*(nIF-1)/2*obj.pitch - offset(1);
                    yIF = linspace(-1,1,nIF)*(nIF-1)/2*obj.pitch - offset(2);
                    [xIF2,yIF2] = ndgrid(xIF,yIF);
                    obj.actuatorCoord = yIF2 + 1i*flip(xIF2);
                    
                    %u0 = ratioTelDm.*linspace(-1,1,resolution)*(nIF-1)/2*(resolution-1)/resolution*obj.pitch;
                    u0 = ratioTelDm.*linspace(-1,1,resolution)*(nIF-1)/2*obj.pitch; % scaled by telescope diamter
                else
                    xIF = 1:nIF;
                    yIF = 1:nIF;
                    [xIF2,yIF2] = ndgrid(xIF,yIF);
                    obj.actuatorCoord = xIF2 + 1i*yIF2;
                    
                    u0 = 1:nIF;
                end
                    
                nValid = sum(validActuator(:));
                kIF = 0;
                
                u = bsxfun( @minus , u0' , xIF );
                wu = zeros(resolution,nIF);
                index_v = u >= -obj.gaussian(end,1) & u <= obj.gaussian(end,1);
                nu = sum(index_v(:));
                wu(index_v) = ppval(obj.splineP,u(index_v));
                
                v = bsxfun( @minus , u0' , yIF);
                wv = zeros(resolution,nIF);
                index_v = v >= -obj.gaussian(end,1) & v <= obj.gaussian(end,1);
                nv = sum(index_v(:));
                wv(index_v) = ppval(obj.splineP,v(index_v));
                
                m_modes = spalloc(resolution^2,nValid,nu*nv);
                
                %             fprintf(' @(influenceFunction)> Computing the 2D DM zonal modes... (%4d,    \n',nValid)
                %             for jIF = 1:nIF
                %
                %                 for iIF = 1:nIF
                %                     if validActuator(iIF,jIF)
                %                         buffer = wv(:,iIF)*wu(:,jIF)';
                %                         kIF = kIF + 1;
                %                         obj.modes(:,kIF) = buffer(:);
                %                                                 fprintf('\b\b\b\b%4d',kIF)
                %                     end
                %                 end
                %
                %             end
                %             fprintf('\n')
                
                indIF = 1:nIF^2;
                indIF(~validActuator) = [];
                [iIF,jIF] = ind2sub([nIF,nIF],indIF);
                kIF = 1:nValid;
                wv = sparse(wv(:,iIF(kIF)));
                wu = sparse(wu(:,jIF(kIF)));
               fprintf(' @(influenceFunction)> Computing the 2D DM zonal modes... (%4d,    \n',nValid)
                for kIF = 1:nValid % parfor doesn't work with sparse matrix!
                    fprintf('\b\b\b\b%4d',kIF)
                    buffer = wv(:,kIF)*wu(:,kIF)';
                    m_modes(:,kIF) = buffer(:);
                end
                fprintf('\n')
                obj.modes = m_modes;
                
                
            else
                
                xIF = real(obj.actuatorCoord(:)') - offset(1);
                yIF = imag(obj.actuatorCoord(:)') - offset(2);
                
                %u0 = linspace(-1,1,resolution)*max(abs(obj.actuatorCoord(:)));
                if ~isempty(diameter)
                    u0 = linspace(-1,1,resolution)*diameter/2; % normalised by the equivalent optical diameter of the DM
                else
                    u0 = linspace(-1,1,resolution)*max(xIF); % normalised by max of the x coordinate of the DM
                end
                nValid = sum(validActuator(:));
                kIF = 0;
                
                u = bsxfun( @minus , u0' , xIF );
                wu = zeros(resolution,nIF);
                index_u = u >= -obj.gaussian(end,1) & u <= obj.gaussian(end,1);
%                 nu = sum(index_u(:));
                wu(index_u) = ppval(obj.splineP,u(index_u));
                
                if ~isempty(diameter)
                    u0 = linspace(-1,1,resolution)*diameter/2; % normalised by the equivalent optical diameter of the DM
                else
                    u0 = linspace(-1,1,resolution)*max(yIF); % normalised by max of the x coordinate of the DM
                end
                v = bsxfun( @minus , u0' , yIF);
                wv = zeros(resolution,nIF);
                index_v = v >= -obj.gaussian(end,1) & v <= obj.gaussian(end,1);
%                 nv = sum(index_v(:));
                wv(index_v) = ppval(obj.splineP,v(index_v));
                expectedNnz = max(sum(index_u))*max(sum(index_v))*nIF;
                add(obj.log,obj,sprintf('Expected non-zeros: %d',expectedNnz))
                add(obj.log,obj,sprintf('Computing the %d 2D DM zonal modes...',nValid))
%                 m_modes = spalloc(resolution^2,nValid,expectedNnz);
                fprintf(' @(influenceFunction)> Computing the 2D DM zonal modes... (%4d,    ',nValid)
                s_i = zeros(expectedNnz,1);
                s_j = zeros(expectedNnz,1);
                s_s = zeros(expectedNnz,1);
                index = 0;
                for kIF = 1:nIF
                    fprintf('\b\b\b\b%4d',kIF)
                    buffer = wv(:,kIF)*wu(:,kIF)';
                    [i_,~,s_] = find(buffer(:));
                    n = length(i_);
                    index = (1:n) + index(end);
                    s_i(index) = i_;
                    s_s(index) = s_;
                    s_j(index) = ones(n,1)*kIF;
%                     m_modes(:,kIF) = buffer(:);
                end
                fprintf('\n')
%                 add(obj.log,obj,sprintf('Actual non-zeros: %d',nnz(m_modes)))
%                 obj.modes = m_modes;
                index = 1:index(end);
                obj.modes = sparse(s_i(index),s_j(index),s_s(index),resolution^2,nValid);
                add(obj.log,obj,sprintf('Actual non-zeros: %d',nnz(obj.modes)))
           end
        end
    end
    
end