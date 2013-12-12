classdef chebtech2 < chebtech
%CHEBTECH2   Approximate smooth functions on [-1,1] with Chebyshev interpolants.
%
%   Class for approximating smooth functions on the interval [-1,1] using
%   function values at 2nd-kind Chebyshev points and coefficients of the
%   corresponding 1st-kind Chebyshev series expansion.
%
% Constructor inputs:
%   CHEBTECH2(OP) constructs a CHEBTECH2 object from the function handle OP. OP
%   should be vectorised (i.e., accept a vector input) and ouput a vector of
%   the same length. CHEBTECH2 objects allow for array-valued construction
%   (i.e., of an array-valued function), in which case OP should accept a column
%   vector of length N and return a matrix of size NxM.
%
%   CHEBTECH2(OP, VSCALE) constructs a CHEBTECH2 with 'happiness' (see
%   HAPPINESSCHECK.m) relative to the maximum of the given vertical scale
%   VSCALE and the (column-wise) infinity norm of the sampled function values
%   of OP, and the fixed horizontal scale HSCALE. If not given (or given as
%   empty), the VSCALE defaults to 0 initially, and HSCALE defaults to 1.
%
%   CHEBTECH2(OP, VSCALE, HSCALE, PREF) overrides the default behavior with
%   that given by the preference structure PREF. See CHEBTECH.PREF for details.
%
%   CHEBTECH2(VALUES, ...) returns a CHEBTECH2 object which interpolates the
%   values in the columns of VALUES at 2nd-kind Chebyshev points and
%   CHEBTECH2({VALUES, COEFFS}, ... ) uses the Chebyshev coefficients passed in
%   COEFFS rather than computing them from VALUES. If VALUES, is empty then the
%   CHEBTECH2 is constructed directly from the COEFFS. If COEFFS are passed,
%   the resulting CHEBTECH2 is always deemed 'happy'.
%
% Examples: % Basic construction: f = chebtech2(@(x) sin(x))
%
%   % Construction with preferences:
%   p.sampleTest = 0; % See CHEBTECH.PREF for details
%   f = chebtech2(@(x) sin(x), [], [], p)
%
%   % Array-valued construction:
%   f = chebtech2(@(x) [sin(x), cos(x), exp(x)])
%
% See also CHEBTECH, CHEBTECH.PREF, CHEBPTS, HAPPINESSCHECK, REFINE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHEBTECH2 Class Description:
%
% The CHEBTECH2 class represents smooth functions on the interval [-1,1] using
% function values at 2nd-kind Chebyshev points and coefficients of the
% corresponding 1st-kind Chebyshev series expansion.
%
% The constructor is supplied with a handle that evaluates a given function on
% an increasingly fine Chebyshev 2nd-kind grid (see REFINE.m) until the
% representation is deemed 'happy' (see HAPPINESSCHECK.m). The resulting object
% can be used to evaluate and operate on the input function.
%
% More information can be found in the CHEBTECH class definition.
%
% Class diagram: [<<CHEBTECH>>] <-- [CHEBTECH2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS IMPLEMENTED BY THIS M-FILE:
    methods
        
        function obj = chebtech2(op, vscale, hscale, pref)
            %Constructor for the CHEBTECH2 class.
            
            % Return an empty CHEBTECH2 on null input:
            if ( (nargin == 0) || isempty(op) )
                return
            end
            
            % Define vscale if none given:
            if ( (nargin < 2) || isempty(vscale) )
                vscale = 0;
            end

            % Define hscale if none given:
            if ( (nargin < 3) || isempty(hscale) )
                hscale = 1;
            end

            % Determine preferences if not given, merge if some are given:
            if ( (nargin < 4) || isempty(pref) )
                pref = chebtech.techPref();
            else
                pref = chebtech.techPref(pref);
            end

            % Force nonadaptive construction if PREF.NUMPOINTS is numeric:
            if ( ~isempty(pref.numPoints) && ~isnan(pref.numPoints) )
                % Evaluate op on the Chebyshev grid of given size:
                op = feval(op, chebtech2.chebpts(pref.numPoints));
            end
            
            % Actual construction takes place here:
            obj = populate(obj, op, vscale, hscale, pref);
            
            if ( obj.ishappy || isnumeric(op) || iscell(op) )
                % No need to error check if we are happy:
                return
            end
            
            % Check for NaNs (if not happy):
            if ( pref.extrapolate )
                % Check for NaNs in interior only (because extrapolate was on):
                if ( any(any(isnan(obj.values(2:end-1,:)))) )
                    error('CHEBFUN:CHEBTECH2:constructor:naneval', ...
                        'Function returned NaN when evaluated.')
                end
                % We make sure not to return NaNs at +1 and -1.
                valuesTemp = extrapolate(obj);
                obj.values([1,end],:) = valuesTemp([1,end],:);
            elseif ( any(isnan(obj.values(:))) )
                % Here we throw an error if NaNs were encountered anywhere.
                error('CHEBFUN:CHEBTECH2:constructor:naneval2', ...
                    'Function returned NaN when evaluated.')
            end
            
        end
        
    end
    
    %% STATIC METHODS IMPLEMENTED BY THIS CLASS:
    methods ( Static = true )
        
        % Aliasing:
        coeffs = alias(coeffs, m)
        
        % Evaluate a Chebyshev interpolant using the barycentric formula:
        out = bary(x, values)
        
        % Compute Chebyshev barycentric weights:
        w = barywts(n)
        
        % Compute Chebyshev points (x) and optionally quadrature (w)
        % and barycentric (v) weights:
        [x, w, v] = chebpts(n);
        
        % Convert coefficients to values:
        values = coeffs2vals(coeffs);
        
        % Make a CHEBTECH2 (constructor shortcut):
        f = make(varargin);
        
        % Compute Chebyshev quadrature weights:
        w = quadwts(n)
        
        % Refinement function for CHEBTECH2 construction (evaluates OP on grid):
        [values, points, giveUp] = refine(op, values, pref)
        
        % Convert values to coefficients:
        coeffs = vals2coeffs(values)
        
    end
    
    %% METHODS IMPLEMENTED BY THIS CLASS:
    methods
        
        % Compose two CHEBTECH2 objects or a CHEBTECH2 with a function handle:
        h = compose(f, op, g, pref)
        
        % Get method:
        val = get(f, prop);
        
    end
    
end
