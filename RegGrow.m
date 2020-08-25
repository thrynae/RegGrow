function [result,OverflowFlag] = RegGrow(I,maxDiff,seed,varargin)
%Region growing algorithm (supports 2D, 3D, and ND)
%
% In each iteration, look at all voxels in the N-neighborhood of the already segmented region. From
% those voxels, add those that are within maxDiff of the mean of the region. Keep repeating the
% loop until no more pixels are added.
% The default neighborhood is the minimal neighborhood (4 for 2D, 6 for 3D, 2*dims for ND)
% As an alternative to entering the neighborhood, a custom kernel can be entered as a binary matrix
% with the same number of dimensions as the input image.
%
% This is slower than a mex implementation would be, but this should be compatible with any release
% and will return the same result on any release (including GNU Octave).
%
%Syntax:
% result = RegGrow()
% result = RegGrow(I)
% result = RegGrow(I,maxDiff)
% result = RegGrow(I,maxDiff,seed)
% result = RegGrow(___,Name,Value)
% result = RegGrow(___,options)
% [result,OverflowFlag]=RegGrow(___)
%
%Note: running without inputs runs on an example image.
%
%output:
% result        The result of the region growing. This will be a logical array of the same size as
%               the input image.
% OverflowFlag  A flag in case of overflow. See below for how to enable checking for an overflow.
%               Normally this output should be false.
%
%input:
%
% I           The image to use, which can be any data type that can be converted to a double.
% maxDiff     The maximum difference between the mean of the currently segmented pixels and
%             candidate pixels that will be added. All candidate pixels are tested at the same
%             time. This value is treated as an absolute value. If omitted, the default value is
%             2/3*std(I(:)).
%             This parameter must be a scalar double, or be convertable to one.
% seed        The initial starting point of the algorithm. If omitted, the default value is the
%             first pixel (linear index 1).
% Name,Value  In addition to the previous parameters, the settings below can be entered with a
%             Name,Value syntax. The parameter names are not case-sensitive.
% options     Instead of the Name,Value, parameters can also be entered in a struct. Missing fields
%             will be set to the default values. Field names are not case-sensitive.
% -nHood      The neighborhood that is checked. By default the minimal neighborhood is used
%             (4-neighborhood for 2D images and 6-neighborhood for 3D images).
%             This setting is ignored if the kernel input is provided.
% -kernel     The kernel input allows more flexibility in the direction of growth. This parameter
%             must be a logical array or be convertible to one. Each dimension must have an odd
%             length. The center pixel is set to true.
%             Note that is parameter describes the direction of growth and is therefore dissimilar
%             to a convolution kernel (i.e. it is flipped along all directions).
% -overflow   This parameter allows the detection of a segmentation leak. This allows the process
%             to exit before it is complete if too many pixels are being segmented. A leak can
%             cause the segmentation to take a very long time to complete, especially for 3D
%             images.
%             If the segmented area or volume exceeds this value (counted as number of pixels or
%             voxels), the second output is set to true and the function halts. No error or warning
%             is triggered.
% -waitbar    Use a waitbar to show the progress of the region growing. The maximum value is either
%             the overflow, or the total number of voxels. This means the waitbar shouldn't reach
%             100%, but does provide a visual indication of progress. The default is false for 2D
%             input and true for higher dimensional input.
% -silent     This is the inverse of the waitbar switch, but it is ignored if wairbar is specified.
%             It is provided to improve code readability.
%
%  _______________________________________________________________________
% | Compatibility | Windows 10  | Ubuntu 20.04 LTS | MacOS 10.15 Catalina |
% |---------------|-------------|------------------|----------------------|
% | ML R2020a     |  works      |  not tested      |  not tested          |
% | ML R2018a     |  works      |  works           |  not tested          |
% | ML R2015a     |  works      |  works           |  not tested          |
% | ML R2011a     |  works      |  works           |  not tested          |
% | ML 6.5 (R13)  |  works      |  not tested      |  not tested          |
% | Octave 5.2.0  |  works      |  works           |  not tested          |
% | Octave 4.4.1  |  works      |  not tested      |  works               |
% """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
%
% Version: 1.1.1
% Date:    2020-08-25
% Author:  H.J. Wisselink
% Licence: CC by-nc-sa 4.0 ( https://creativecommons.org/licenses/by-nc-sa/4.0 )
% Email = 'h_j_wisselink*alumnus_utwente_nl';
% Real_email = regexprep(Email,{'*','_'},{'@','.'})

% Setting the shift array is also possible, but is undocumented and will not be checked when
% parsing the parameters. The shift array will be converted to a kernel.
if nargout>2
    error('HJW:RegGrow:nargout','Incorrect number of output arguments.')
end

inputs=cell(1,nargin);
if nargin>=1,inputs{1}=I;end
if nargin>=2,inputs{2}=maxDiff;end
if nargin>=3,inputs{3}=seed;end
if nargin>=4,inputs(4:end)=varargin;end
[success,parsed,ME]=RegGrow_ParseInputs(inputs{:});
if ~success
    % The throwAsCaller function was introduced in R2007b, hence the rethrow here.
    rethrow(ME)
else
    [I,maxDiff,seed,SE,overflow,wait]=deal(...
        parsed.I,parsed.maxDiff,parsed.seed,parsed.kernel,parsed.overflow,parsed.waitbar);
end

if wait.bar
    wait.h=waitbar(0,'processing...');
    wait.count=0;%initialize
	n_extra_dots=0;
end

OverflowFlag=false;%initialize overflow parameter
result=logical(zeros(size(I))); %#ok<LOGL>
seed=num2cell(seed);result(seed{:})=true;%pre-allocate the result
newpoints=true;%initialize to start loop
while newpoints
    if overflow.check
        if sum(result(:))>overflow.volume
            OverflowFlag=true;
            break
        end
    end
    
    if wait.bar
        wait.count=wait.count+1;
        if mod(wait.count,5)==0
            try
                n_extra_dots=mod(n_extra_dots+1,3);
                dots=repmat('.',1,1+n_extra_dots);
                txt=sprintf('processing%s',dots);
                waitbar(sum(result(:))/wait.max,wait.h,txt),drawnow
            catch
                ME=lasterror; %#ok<LERR>
                if strcmp(ME.identifier,'MATLAB:waitbar:InvalidSecondInput')
                    % The waitbar was closed, so now we have the choice of interpreting this as
                    % either a forced quit, or a mistake. In the first case, an error should be
                    % thrown, in the second case wait.bar is set to false.
                    wait.bar=false;
                end
            end
        end
    end
    
    % Get the shell by using a convolution (whose result is a double array of zeros and ones), then
    % mark all voxels that have a value outside intensity limits as false. This has to be done
    % every loop, as the intensity limits can shift when adding voxels. This syntax is equivalent
    % to or(candidates,reject). Calculating the full reject matrix is probably more efficient than
    % calculating it only for the candidates and then having to index into it.
    
    candidates= ~result & convn(result,SE,'same')>0.5;%shell
    reject=abs(I-mean(I(result)))>=maxDiff;
    candidates(reject)=false;
    if any(candidates(:))
        result(candidates)=true;
    else
        newpoints=false;
    end
end
if wait.bar
    delete(wait.h)
end
if nargout==0
    %Prevent cluttering of the command window / workspace.
    clear result
end
end
function [success,options,ME]=RegGrow_ParseInputs(varargin)
%Parse the inputs. See top function for syntax options and defaults.
%Input of the shift array is undocumented (and is overwritten by both the nHood and kernel
%parameters).

%pre-assign outputs
success=false;
options=struct;
ME=struct('identifier','','message','');

persistent default
if isempty(default)
    % The default example image is taken from the image generated by the image function. For the
    % full back story, see this blog post:
    %     https://blogs.mathworks.com/steve/2006/10/17/the-story-behind-the-matlab-default-image/
    %
    % In GNU Octave this image is not available and it is replaced by a logo image. The example
    % setup will still work as intended.
    isOctave=exist('OCTAVE_VERSION', 'builtin');
    if isOctave
        IM = get(0,'DefaultImageCData');
        IM = IM/max(IM(:));
    else
        % Because the behavior of bitcmp changed in R2014b, the bitslice function included in the
        % blog post no longer works. The code below should work on all releases again.
        defimage = pow2(get(0,'DefaultImageCData'),47);
        IM = bitshift(defimage,-37);
        IM = fix(IM);
        IM = bitand(IM,31);%bitcmp(0,numbits) evaluates to 31 here
        IM = IM/max(IM(:));
    end
    
    default.I=IM;
    default.maxDiff=[];%dynamic default: std(I(:))
    default.seed=[];%dynamic default: first position
    default.nHood=[];%4 neighborhood for 2D, 6 for 3D
    default.kernel=[];
    default.shifts=[];%dynamic default: based on image and nHood, or kernel
    default.overflow=[];%dynamic default: struct with .check and .volume
    default.KeepShellVoxelsOnly=false;
    default.waitbar=[];%dynamic default: struct with .bar (false for 2D)
end

if nargin==0
    options=RegGrow_check_dynamic_defaults(default);
    options.seed=[2 2];%the example image on Matlab has a border
    success=true;
    return
end

%Parse the first inputs (I, maxDiff, and seed).
counter=0;names={'I','maxDiff','seed'};
while numel(varargin)>0 && ~(isa(varargin{1},'char') || isa(varargin{1},'struct'))
    counter=counter+1;
    options.(names{counter})=varargin{1};
    %Remove the current entry for next loop, or for parsing of the options struct or the Name,Value
    %pairs.
    varargin(1)=[];
end

if numel(varargin)>0
    %Now only the Name,Value pairs or the options struct remain.
    if isa(varargin{1},'struct')
        %I, maxDiff, and seed might already be loaded
        opts=options;
        options=varargin{1};
        %merge the two structs
        fn=fieldnames(opts);
        for n=1:numel(fn)
            options.(fn{n})=opts.(fn{n});
        end
    else
        try
            %Loop through the Name,Value pairs.
            while numel(varargin)>0
                options.(varargin{1})=varargin{2};
                varargin(1:2)=[];
            end
        catch
            ME.identifier='HJW:RegGrow:ParseError';
            ME.message='Parsing of the inputs failed.';
            return
        end
    end
end

%The inputs are stored in the struct, so the inputs can be checked. The check is not
%case-sensitive. To ensure proper capitalization the parameters are explicitly stored.
fn=fieldnames(options);
[ignore,order]=sort(lower(fn));fn=fn(order); %#ok<ASGLU>
for k=1:numel(fn)
    curr_option=fn{k};
    item=options.(curr_option);
    ME.identifier=['HJW:RegGrow:incorrect_input_opt_' lower(curr_option)];
    switch lower(curr_option)
        %NB: If a new parameter is alphabetically before the image parameter, it can't rely on IM
        %already being a double, so special attention is necesarry for user-defined classes.
        case 'i'
            try
                item=double(item);
            catch
                ME.message='The image input must be convertable to a double .';
                return
            end
            options.I=item;
        case 'kernel'
            try
                %If the image provided as an input, it is already loaded to the the options struct.
                if isfield(options,'I'),IM=options.I;else,IM=default.I;end
                item=logical(item);
                if ndims(item)~=ndims(IM) || ~all(mod(size(item),2)==1)
                    disp(item(-1))%trigger error
                end
                tmp=item;center=(size(tmp)+1)/2;center=num2cell(center);tmp(center{:})=false;
                item(center{:})=true;%doesn't really affect calculation
                if sum(tmp(:))==0
                    warning('HJW:RegGrow:CentroidOnlyKernel',...
                        ['Only the centroid of the kernel is marked as true.\n',...
                        'This will cause the result to only contain the seed position.'])
                    %Technically this message will trigger regardless of the actual value of the
                    %centroid. Because the value of the centroid is always set to true, this is a
                    %distinction without a difference. To keep the message easy to read that point
                    %is omitted.
                end
            catch
                ME.message=['The kernel input must a logical with ndims(IM) dimensions.',...
                    char(10),'Each dimension must be an odd length.']; %#ok<CHARTEN>
                return
            end
            options.kernel=item;
        case 'maxdiff'
            try item=double(item);catch,end %an error here will trigger an error on the next line
            if ~isa(item,'double') || numel(item)~=1 || item<0 || isnan(item)
                ME.message=['The maxDiff input must be a non-negative numeric scalar.',char(10),...
                    'It must also be convertable to a double.']; %#ok<CHARTEN>
                return
            end
            options.maxDiff=item;
        case 'nhood'
            %Skip if kernel or shifts is set
            if any(ismember(lower(fn),{'kernel','shifts'}))
                continue
            end
            %If the image provided as an input, it is already loaded to the the options struct.
            if isfield(options,'I'),IM=options.I;else,IM=default.I;end
            
            msg=['The nHood parameter must be either the maximal or minimal value, ',...
                char(10),'so 4 or 8 for 2D images, 6 or 26 for 3D images, etc.',char(10),...
                'This parameter must be set if the kernel parameter is not set.']; %#ok<CHARTEN>
            %These tests will not trigger an error if this parameter is set to an empty array.
            try
                dims=ndims(IM);
                if ~ismember(item,[3^dims-1 2*dims])%maximal and minimal
                    ME.message=msg;
                    return
                end
            catch
                ME.message=msg;
                return
            end
            options.nHood=item;
        case 'overflow'
            if ~isempty(item) && ( ~isnumeric(item) || numel(item)>1 || abs(item-round(item))>eps)
                ME.message='The overflow volume must be empty, or a numeric scalar integer.';
                return
            end
            options.overflow=item;
        case 'seed'
            %If the image provided as an input, it is already loaded to the the options struct.
            if isfield(options,'I'),IM=options.I;else,IM=default.I;end
            try
                if numel(item)~=ndims(IM)
                    ME.message='The seed must be a vector with a valid position.';
                    return
                end
                seed=item(:)';
                seed_=num2cell(seed);
                val=IM(seed_{:});%#ok<NASGU>
                %if this works, than the seed must be a valid vector
            catch
                return
            end
            ME.message='';
            options.seed=seed;
        case {'waitbar','silent'} %same parameter, two names
            [passed,item]=test_if_scalar_logical(item);
            if ~passed
                ME.message='waitbar must be a logical scalar.';
                return
            end
            %invert if 'silent' was entered
            if strcmpi(curr_option,'silent'),item= ~item;end
            options.waitbar=item;
        case 'shifts'
            options.shifts=item;%undocumented, don't do any checks
        otherwise
            ME.message=sprintf('Name,Value pair not recognized: %s',curr_option);
            ME.identifier='HJW:RegGrow:incorrect_input_NameValue';
            return
    end
end

%Now the explicit inputs are loaded, so the defaults can be added for any missing parameters.
fn=fieldnames(default);
for k=1:numel(fn)
    if ~isfield(options,fn(k))
        options.(fn{k})=default.(fn{k});
    end
end

%Set the success flag and check the dynamic defaults before returning.
success=true;
options=RegGrow_check_dynamic_defaults(options);
end
function options=RegGrow_check_dynamic_defaults(options)
%set the dynamic defaults if the options were not set.
if isempty(options.seed)
    %create a vector of ones
    s=ones(ndims(options.I),1);
    options.seed=s;
end
if isempty(options.maxDiff)
    options.maxDiff=2/3*std(options.I(:));
end

if ~isempty(options.shifts)%override both kernel and nHood inputs
    %construct the kernel from the shifts variable
    shifts=options.shifts;
    sz=max(abs(shifts(:)));
    SE=zeros( (1+2*sz)*ones(1,size(shifts,2)) );
    shifts=shifts+sz+1;
    shifts=mat2cell(shifts,size(shifts,1),ones(1,size(shifts,2)));
    ind=sub2ind(size(SE),shifts{:});
    SE(ind)=1;SE=logical(SE);
    options.kernel=SE;
end
if isempty(options.kernel)
    nHood=options.nHood;
    %set default if missing (minimal connectivity)
    if isempty(nHood)
        nHood=2*ndims(options.I);
    end
    dims=ndims(options.I);
    if nHood==(2*dims) %minimal
        %create a sphere with r=1
        XYZ=cell(dims,1);
        [XYZ{:}]=ndgrid(-1:1);
        XYZ=cat(dims+1,XYZ{:});
        XYZ=sqrt(sum(XYZ.^2,dims+1));
        SE=XYZ<=1;
    elseif nHood==(3^dims-1) %maximal
        SE=logical(ones(3*ones(1,dims))); %#ok<LOGL>
        %the doc for conndef describes this as ones(repmat(3,1,NUM_DIMS))
    end
    options.kernel=SE;
end
val=options.overflow;
if isempty(val) || isnan(val) || isinf(val)
    options.overflow=struct('check',false,'volume',inf);
else
    options.overflow=struct('check',true,'volume',val);
end
if isempty(options.waitbar)
    options.waitbar= ndims(options.I)~=2 ; %#ok<ISMAT>
end
options.waitbar=struct('bar',options.waitbar,'max',min(numel(options.I),options.overflow.volume));
end
function [passed,item]=test_if_scalar_logical(item)
%test if the input is a scalar logical or convertable to it
%(use isLogical to trigger an input error, use val as the parsed input)
%
% Allowed values:
%- true or false
%- 1 or 0
%- 'on' or 'off'
persistent states
if isempty(states)
    states={true,false;...
        1,0;...
        'on','off'};
    try
        states(end+1,:)=eval('{"on","off"}');
    catch
    end
end
passed=true;
try
    for n=1:size(states,1)
        for m=1:2
            if isequal(item,states{n,m})
                item=states{1,m};return
            end
        end
    end
    if isa(item,'matlab.lang.OnOffSwitchState')
        item=logical(item);return
    end
catch
end
passed=false;
end