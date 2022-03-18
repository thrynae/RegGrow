function [result,OverflowFlag] = RegGrow(varargin)
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
% Syntax:
%   result = RegGrow()
%   result = RegGrow(I)
%   result = RegGrow(I,maxDiff)
%   result = RegGrow(I,maxDiff,seed)
%   result = RegGrow(___,Name,Value)
%   result = RegGrow(___,options)
%   [result,OverflowFlag]=RegGrow(___)
%Note: running without inputs runs on an example image.
%
% Input/output arguments:
% result:
%   The result of the region growing. This will be a logical array of the same size as the input
%   image.
% OverflowFlag:
%   A flag in case of overflow. See below for how to enable checking for an overflow. Normally this
%   output should be false.
% I:
%   The image to use, which can be any data type that can be converted to a double.
% maxDiff:
%   The maximum difference between the mean of the currently segmented pixels and candidate pixels
%   that will be added. All candidate pixels are tested at the same time. This value is treated as
%   an absolute value. If omitted, the default value is 2/3*std(I(:)).
%   This parameter must be a scalar double, or be convertible to one.
% seed:
%   The initial starting point of the algorithm. If omitted, the default value is the first pixel
%   (linear index 1).
% options:
%   A struct with Name,Value parameters. Missing parameters are filled with the defaults listed
%   below. Using incomplete parameter names or incorrect capitalization is allowed, as long as
%   there is a unique match.
%
% Name,Value parameters:
%   nHood:
%      The neighborhood that is checked. By default the minimal neighborhood is used
%      (4-neighborhood for 2D images and 6-neighborhood for 3D images).
%      This setting is ignored if the kernel input is provided.
%   kernel:
%      The kernel input allows more flexibility in the direction of growth. This parameter must be
%      a logical array or be convertible to one. Each dimension must have an odd length. The center
%      pixel is set to true. Note that is parameter describes the direction of growth and is
%      therefore dissimilar to a convolution kernel (i.e. it is flipped along all directions).
%   overflow:
%      This parameter allows the detection of a segmentation leak. This allows the process to exit
%      before it is complete if too many pixels are being segmented. A leak can cause the
%      segmentation to take a very long time to complete, especially for 3D images.
%      If the segmented area or volume exceeds this value (counted as number of pixels or voxels),
%      the second output is set to true and the function halts. No error or warning is triggered.
%   waitbar:
%      Use a waitbar to show the progress of the region growing. The maximum value is either the
%      overflow, or the total number of voxels. This means the waitbar shouldn't reach 100%, but
%      does provide a visual indication of progress. The default is false for 2D input and true for
%      higher dimensional input.
%   silent:
%      This is the inverse of the waitbar switch, but it is ignored if wairbar is specified. It is
%      provided to improve code readability.
%   WaitbarClosingBehavior:
%      This parameter must be 'warning', 'error', 'return', or 'ignore'. The default is 'ignore'.
%      If set to warn or error, closing the waitbar will cause a warning/error to be triggered.
%      Setting to 'return' will cause the function to return the incomplete result.
%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%|                                                                         |%
%|  Version: 1.2.0                                                         |%
%|  Date:    2022-03-18                                                    |%
%|  Author:  H.J. Wisselink                                                |%
%|  Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 ) |%
%|  Email = 'h_j_wisselink*alumnus_utwente_nl';                            |%
%|  Real_email = regexprep(Email,{'*','_'},{'@','.'})                      |%
%|                                                                         |%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%
% Tested on several versions of Matlab (ML 6.5 and onward) and Octave (4.4.1 and onward), and on
% multiple operating systems (Windows/Ubuntu/MacOS). For the full test matrix, see the HTML doc.
% Compatibility considerations:
% - This is expected to work on all releases.


% Setting the shift array is also possible, but is undocumented and will not be checked when
% parsing the parameters. The shift array will be converted to a kernel.
if nargout>2
    error('HJW:RegGrow:nargout','Incorrect number of output arguments.')
end

[success,parsed,ME]=RegGrow_ParseInputs(varargin{:});
if ~success
    % The throwAsCaller function was introduced in R2007b, hence the rethrow here.
    rethrow(ME)
else
    [I,maxDiff,seed,SE,overflow,wait]=deal(...
        parsed.I,parsed.maxDiff,parsed.seed,parsed.kernel,parsed.overflow,parsed.waitbar);
end

if wait.bar
    wait.h=waitbar(0,'processing...');
    wait.count=0;%Initialize.
    n_extra_dots=0;
end

OverflowFlag=false; % Initialize overflow parameter.
result=logical(zeros(size(I))); %#ok<LOGL>
seed=num2cell(seed);result(seed{:})=true; % Pre-allocate the result.
newpoints=true; % Initialize to start loop.
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
                    % The waitbar was closed.
                    switch parsed.WaitbarClosingBehavior
                        case 'warning'
                            warning('HJW:RegGrow:WaitbarClosedWarning',...
                                'The waitbar was closed, continuing...')
                            wait.bar=false;
                        case 'error'
                            error('HJW:RegGrow:WaitbarClosedError',...
                                'The waitbar was closed.')
                        case 'return'
                            if nargout==0
                                %Prevent cluttering of the command window / workspace.
                                clear(var2str(result))
                            end
                            return
                        case 'ignore'
                            wait.bar=false;
                    end
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
    clear(var2str(result))
end
end
function tf=ifversion(test,Rxxxxab,Oct_flag,Oct_test,Oct_ver)
%Determine if the current version satisfies a version restriction
%
% To keep the function fast, no input checking is done. This function returns a NaN if a release
% name is used that is not in the dictionary.
%
% Syntax:
%   tf=ifversion(test,Rxxxxab)
%   tf=ifversion(test,Rxxxxab,'Octave',test_for_Octave,v_Octave)
%
% Input/output arguments:
% tf:
%   If the current version satisfies the test this returns true. This works similar to verLessThan.
% Rxxxxab:
%   A char array containing a release description (e.g. 'R13', 'R14SP2' or 'R2019a') or the numeric
%   version (e.g. 6.5, 7, or 9.6).
% test:
%   A char array containing a logical test. The interpretation of this is equivalent to
%   eval([current test Rxxxxab]). For examples, see below.
%
% Examples:
% ifversion('>=','R2009a') returns true when run on R2009a or later
% ifversion('<','R2016a') returns true when run on R2015b or older
% ifversion('==','R2018a') returns true only when run on R2018a
% ifversion('==',9.9) returns true only when run on R2020b
% ifversion('<',0,'Octave','>',0) returns true only on Octave
% ifversion('<',0,'Octave','>=',6) returns true only on Octave 6 and higher
%
% The conversion is based on a manual list and therefore needs to be updated manually, so it might
% not be complete. Although it should be possible to load the list from Wikipedia, this is not
% implemented.
%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%|                                                                         |%
%|  Version: 1.1.1                                                         |%
%|  Date:    2022-03-14                                                    |%
%|  Author:  H.J. Wisselink                                                |%
%|  Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 ) |%
%|  Email = 'h_j_wisselink*alumnus_utwente_nl';                            |%
%|  Real_email = regexprep(Email,{'*','_'},{'@','.'})                      |%
%|                                                                         |%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%
% Tested on several versions of Matlab (ML 6.5 and onward) and Octave (4.4.1 and onward), and on
% multiple operating systems (Windows/Ubuntu/MacOS). For the full test matrix, see the HTML doc.
% Compatibility considerations:
% - This is expected to work on all releases.

% The decimal of the version numbers are padded with a 0 to make sure v7.10 is larger than v7.9.
% This does mean that any numeric version input needs to be adapted. multiply by 100 and round to
% remove the potential for float rounding errors.
% Store in persistent for fast recall (don't use getpref, as that is slower than generating the
% variables and makes updating this function harder).
persistent  v_num v_dict octave
if isempty(v_num)
    % Test if Octave is used instead of Matlab.
    octave=exist('OCTAVE_VERSION', 'builtin');
    
    % Get current version number. This code was suggested by Jan on this thread:
    % https://mathworks.com/matlabcentral/answers/1671199#comment_2040389
    v_num = [100, 1] * sscanf(version, '%d.%d', 2);
    
    % Get dictionary to use for ismember.
    v_dict={...
        'R13' 605;'R13SP1' 605;'R13SP2' 605;'R14' 700;'R14SP1' 700;'R14SP2' 700;
        'R14SP3' 701;'R2006a' 702;'R2006b' 703;'R2007a' 704;'R2007b' 705;
        'R2008a' 706;'R2008b' 707;'R2009a' 708;'R2009b' 709;'R2010a' 710;
        'R2010b' 711;'R2011a' 712;'R2011b' 713;'R2012a' 714;'R2012b' 800;
        'R2013a' 801;'R2013b' 802;'R2014a' 803;'R2014b' 804;'R2015a' 805;
        'R2015b' 806;'R2016a' 900;'R2016b' 901;'R2017a' 902;'R2017b' 903;
        'R2018a' 904;'R2018b' 905;'R2019a' 906;'R2019b' 907;'R2020a' 908;
        'R2020b' 909;'R2021a' 910;'R2021b' 911;'R2022a' 912};
end

if octave
    if nargin==2
        warning('HJW:ifversion:NoOctaveTest',...
            ['No version test for Octave was provided.',char(10),...
            'This function might return an unexpected outcome.']) %#ok<CHARTEN>
        if isnumeric(Rxxxxab)
            v=0.1*Rxxxxab+0.9*fix(Rxxxxab);v=round(100*v);
        else
            L=ismember(v_dict(:,1),Rxxxxab);
            if sum(L)~=1
                warning('HJW:ifversion:NotInDict',...
                    'The requested version is not in the hard-coded list.')
                tf=NaN;return
            else
                v=v_dict{L,2};
            end
        end
    elseif nargin==4
        % Undocumented shorthand syntax: skip the 'Octave' argument.
        [test,v]=deal(Oct_flag,Oct_test);
        % Convert 4.1 to 401.
        v=0.1*v+0.9*fix(v);v=round(100*v);
    else
        [test,v]=deal(Oct_test,Oct_ver);
        % Convert 4.1 to 401.
        v=0.1*v+0.9*fix(v);v=round(100*v);
    end
else
    % Convert R notation to numeric and convert 9.1 to 901.
    if isnumeric(Rxxxxab)
        v=0.1*Rxxxxab+0.9*fix(Rxxxxab);v=round(100*v);
    else
        L=ismember(v_dict(:,1),Rxxxxab);
        if sum(L)~=1
            warning('HJW:ifversion:NotInDict',...
                'The requested version is not in the hard-coded list.')
            tf=NaN;return
        else
            v=v_dict{L,2};
        end
    end
end
switch test
    case '==', tf= v_num == v;
    case '<' , tf= v_num <  v;
    case '<=', tf= v_num <= v;
    case '>' , tf= v_num >  v;
    case '>=', tf= v_num >= v;
end
end
function [opts,replaced]=parse_NameValue(default,varargin)
%Match the Name,Value pairs to the default option, ignoring incomplete names, case, and underscores
%
% The first output is a struct with the same fields as the first input, with field contents
% replaced according to the supplied options struct or Name,Value pairs.
% The second output is a cellstr containing the field names that have been set.
%
% If this fails to find a match, this will throw an error with the offending name as the message.
%
% If there are multiple occurences of a Name, only the last Value will be returned. This is the
% same as Matlab internal functions like plot. GNU Octave also has this behavior.
%
% If a struct array is provided, only the first element will be used. An empty struct array will
% trigger an error.

switch numel(default)
    case 0
        error('parse_NameValue:MixedOrBadSyntax',...
            'Optional inputs must be entered as Name,Value pairs or as a scalar struct.')
    case 1
        % Do nothing.
    otherwise
        % If this is a struct array, explicitly select the first element.
        default=default(1);
end

% Create default output and return if no other inputs exist.
opts=default;replaced={};
if nargin==1,return,end

% Unwind an input struct to Name,Value pairs.
try
    struct_input=numel(varargin)==1 && isa(varargin{1},'struct');
    NameValue_input=mod(numel(varargin),2)==0 && all(...
        cellfun('isclass',varargin(1:2:end),'char'  ) | ...
        cellfun('isclass',varargin(1:2:end),'string')   );
    if ~( struct_input || NameValue_input )
        error('trigger')
    end
    if nargin==2
        Names=fieldnames(varargin{1});
        Values=struct2cell(varargin{1});
    else
        % Wrap in cellstr to account for strings (this also deals with the fun(Name=Value) syntax).
        Names=cellstr(varargin(1:2:end));
        Values=varargin(2:2:end);
    end
    if ~iscellstr(Names),error('trigger');end %#ok<ISCLSTR>
catch
    % If this block errors, that is either because a missing Value with the Name,Value syntax, or
    % because the struct input is not a struct, or because an attempt was made to mix the two
    % styles. In future versions of this functions an effort might be made to handle such cases.
    error('parse_NameValue:MixedOrBadSyntax',...
        'Optional inputs must be entered as Name,Value pairs or as a scalar struct.')
end

% The fieldnames will be converted to char matrices in the section below. First an exact match is
% tried, then a case-sensitive (partial) match, then ignoring case, followed by ignoring any
% underscores.
default_Names=fieldnames(default);
Names_char=cell(1,3);
Names_cell{1}=default_Names;
Names_cell{2}=lower(Names_cell{1});
Names_cell{3}=strrep(Names_cell{2},'_','');

% Attempt to match the names.
replaced=false(size(default_Names));
for n=1:numel(Names)
    name=Names{n};
    
    % Try a case-sensitive match.
    [match_idx,Names_char{1}]=parse_NameValue__find_match(Names_char{1},Names_cell{1},name);
    
    % Try a case-insensitive match.
    if numel(match_idx)~=1
        name=lower(name);
        [match_idx,Names_char{2}]=parse_NameValue__find_match(Names_char{2},Names_cell{2},name);
    end
    
    % Try a case-insensitive match ignoring underscores.
    if numel(match_idx)~=1
        name=strrep(name,'_','');
        [match_idx,Names_char{3}]=parse_NameValue__find_match(Names_char{3},Names_cell{3},name);
    end
    
    if numel(match_idx)~=1
        error('parse_NameValue:NonUniqueMatch',Names{n})
    end
    
    %Store the Value in the output struct and mark it as replaced.
    opts.(default_Names{match_idx})=Values{n};
    replaced(match_idx)=true;
end
replaced=default_Names(replaced);
end
function [match_idx,Names_char]=parse_NameValue__find_match(Names_char,Names_cell,name)
% Try to match the input field to the fields of the struct.

% First attempt an exact match.
match_idx=find(ismember(Names_cell,name));
if numel(match_idx)==1,return,end

% Only spend time building the char array if this point is reached.
if isempty(Names_char),Names_char=parse_NameValue__name2char(Names_cell);end

% Since the exact match did not return a unique match, attempt to match the start of each array.
% Select the first part of the array. Since Names is provided by the user it might be too long.
tmp=Names_char(:,1:min(end,numel(name)));

% Find the number of non-matching characters on every row. The cumprod on the logical array is
% to make sure that only the starting match is considered.
non_matching=numel(name)-sum(cumprod(double(tmp==repmat(name,size(tmp,1),1)),2),2);
match_idx=find(non_matching==0);
end
function Names_char=parse_NameValue__name2char(Names_char)
% Convert a cellstr to a padded char matrix.
len=cellfun('prodofsize',Names_char);maxlen=max(len);
for n=find(len<maxlen).' % Pad with spaces where needed
    Names_char{n}((end+1):maxlen)=' ';
end
Names_char=vertcat(Names_char{:});
end
function [success,options,ME]=RegGrow_ParseInputs(varargin)
% Parse the inputs. See top function for syntax options and defaults.
% Input of the shift array is undocumented (and is overwritten by both the nHood and kernel
% parameters).

% Pre-assign outputs
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
    if ifversion('<',0,'Octave','>',0)
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
    default.silent=[];%merged into waitbar
    default.ErrorOnWaitbarClosing=false;
end

if nargin==0
    options=RegGrow_check_dynamic_defaults(default);
    options.seed=[2 2]; % The example image on Matlab has a border.
    success=true;
    return
end

% Parse the first inputs (I, maxDiff, and seed).
counter=0;names={'I','maxDiff','seed'};
while numel(varargin)>0 && ~(isa(varargin{1},'char') || isa(varargin{1},'struct'))
    counter=counter+1;
    options.(names{counter})=varargin{1};
    % Remove the current entry for next loop, or for parsing of the options struct or the
    % Name,Value pairs.
    varargin(1)=[];
end
try ME_=[];[options,replaced]=parse_NameValue(default,options); %#ok<NASGU>
catch ME_;if isempty(ME_),ME_=lasterror;end,ME=ME_;return,end %#ok<LERR>

% Now only the Name,Value pairs or the options struct remain.
if numel(varargin)>0
    try ME_=[];[options,replaced2]=parse_NameValue(options,varargin{:}); %#ok<NASGU>
    catch ME_;if isempty(ME_),ME_=lasterror;end,ME=ME_;return,end %#ok<LERR>
    
    % The Name,Value pairs or option struct might specify I, maxDiff, or seed.
    replaced=unique([replaced;replaced2]);
end

% Loop through the replaced options in alphabetical order.
[ignore,order]=sort(lower(replaced));replaced=replaced(order); %#ok<ASGLU>
for k=1:numel(replaced)
    curr_option=replaced{k};
    item=options.(curr_option);
    ME.identifier=['HJW:RegGrow:incorrect_input_opt_' lower(curr_option)];
    switch lower(curr_option)
        %NB: If a new parameter is alphabetically before the image parameter, it can't rely on IM
        %already being a double, so special attention is necessary for user-defined classes.
        case 'i'
            try
                item=double(item);
            catch
                ME.message='The image input must be convertible to a double.';
                return
            end
            options.I=item;
        case 'kernel'
            try
                % If the image provided as an input, it is already loaded to the options struct.
                IM=options.I;
                item=logical(item);
                if ndims(item)~=ndims(IM) || ~all(mod(size(item),2)==1)
                    disp(item(-1)) % Trigger error.
                end
                tmp=item;center=(size(tmp)+1)/2;center=num2cell(center);tmp(center{:})=false;
                item(center{:})=true; % The center doesn't really affect calculation.
                if sum(tmp(:))==0
                    warning('HJW:RegGrow:CentroidOnlyKernel',...
                        ['Only the centroid of the kernel is marked as true.\n',...
                        'This will cause the result to only contain the seed position.'])
                    % Technically this message will trigger regardless of the actual value of the
                    % centroid. Because the value of the centroid is always set to true, this is a
                    % distinction without a difference. To keep the message easy to read that point
                    % is omitted.
                end
            catch
                ME.message=['The kernel input must a logical with ndims(IM) dimensions.',...
                    char(10),'Each dimension must be an odd length.']; %#ok<CHARTEN>
                return
            end
            options.kernel=item;
        case 'maxdiff'
            try item=double(item);catch,end % An error here will trigger an error on the next line.
            if ~isa(item,'double') || numel(item)~=1 || item<0 || isnan(item)
                ME.message=['The maxDiff input must be a non-negative numeric scalar.',char(10),...
                    'It must also be convertible to a double.']; %#ok<CHARTEN>
                return
            end
            options.maxDiff=item;
        case 'nhood'
            % Skip if kernel or shifts is set.
            if any(ismember(lower(replaced),{'kernel','shifts'}))
                continue
            end
            % If the image provided as an input, it is already loaded to the options struct.
            IM=options.I;
            
            msg=['The nHood parameter must be either the maximal or minimal value, ',...
                char(10),'so 4 or 8 for 2D images, 6 or 26 for 3D images, etc.',char(10),...
                'This parameter must be set if the kernel parameter is not set.']; %#ok<CHARTEN>
            % These tests will not trigger an error if this parameter is set to an empty array.
            try
                dims=ndims(IM);
                if numel(item)>1,item=0;end
                item=double(item);
                if ~ismember(item,[3^dims-1 2*dims]) % Maximal and minimal.
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
            % If the image provided as an input, it is already loaded to the options struct.
            IM=options.I;
            msg='The seed must be a vector with a valid position.';
            try
                if numel(item)~=ndims(IM)
                    ME.message=msg;
                    return
                end
                seed=item(:)';
                seed_=num2cell(seed);
                val=IM(seed_{:});%#ok<NASGU>
                % If this worked the seed must be a valid vector.
            catch
                ME.message=msg;
                return
            end
            ME.message='';
            options.seed=seed;
        case {'waitbar','silent'} % Same parameter, two names.
            [passed,item]=test_if_scalar_logical(item);
            if ~passed
                ME.message='waitbar must be a logical scalar.';
                return
            end
            % Invert if 'silent' was entered.
            if strcmpi(curr_option,'silent'),item= ~item;end
            options.waitbar=item;
        case 'shifts'
            options.shifts=item; % Undocumented, don't do any checks.
        case 'waitbarclosingbehavior'
            s=struct('warning',0,'error',0,'return',0,'ignore',0);
            [ignore,item]=parse_NameValue(s,{item,true}); %#ok<ASGLU> 
            if ~passed
                ME.message=['WaitbarClosingBehavior must be one of these:' char(10),...
                    '''warning'', ''error'', ''return'', or ''ignore''.']; %#ok<CHARTEN> 
                return
            end
            options.WaitbarClosingBehavior=item;
    end
end

% Set the success flag and check the dynamic defaults before returning.
success=true;
options=RegGrow_check_dynamic_defaults(options);
end
function options=RegGrow_check_dynamic_defaults(options)
% Set the dynamic defaults if the options were not set.
if isempty(options.seed)
    % Create a vector of ones.
    s=ones(ndims(options.I),1);
    options.seed=s;
end
if isempty(options.maxDiff)
    options.maxDiff=2/3*std(options.I(:));
end

if ~isempty(options.shifts) % Override both kernel and nHood inputs.
    % Construct the kernel from the shifts variable.
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
    % Set default if missing (minimal connectivity).
    if isempty(nHood)
        nHood=2*ndims(options.I);
    end
    dims=ndims(options.I);
    if nHood==(2*dims) %Minimal
        % Create a sphere with r=1.
        XYZ=cell(dims,1);
        [XYZ{:}]=ndgrid(-1:1);
        XYZ=cat(dims+1,XYZ{:});
        XYZ=sqrt(sum(XYZ.^2,dims+1));
        SE=XYZ<=1;
    elseif nHood==(3^dims-1) %Maximal
        SE=logical(ones(3*ones(1,dims))); %#ok<LOGL>
        % The doc for conndef describes this as ones(repmat(3,1,NUM_DIMS)).
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
    options.waitbar= ndims(options.I)>2 ; %#ok<ISMAT>
end
options.waitbar=struct('bar',options.waitbar,'max',min(numel(options.I),options.overflow.volume));
end
function [isLogical,val]=test_if_scalar_logical(val)
%Test if the input is a scalar logical or convertible to it.
%The char and string test are not case sensitive.
%(use the first output to trigger an input error, use the second as the parsed input)
%
% Allowed values:
%- true or false
%- 1 or 0
%- 'on' or 'off'
%- "on" or "off"
%- matlab.lang.OnOffSwitchState.on or matlab.lang.OnOffSwitchState.off
%- 'enable' or 'disable'
%- 'enabled' or 'disabled'
persistent states
if isempty(states)
    states={true,false;...
        1,0;...
        'on','off';...
        'enable','disable';...
        'enabled','disabled'};
    if ifversion('>=','R2016b','<',0)
        states(end+1,:)={string('on'),string('off')}; %#ok<STRQUOT>
    end
end
isLogical=true;
try
    if isa(val,'char') || isa(val,'string')
        try val=lower(val);catch,end
    end
    for n=1:size(states,1)
        for m=1:2
            if isequal(val,states{n,m})
                val=states{1,m};return
            end
        end
    end
    if isa(val,'matlab.lang.OnOffSwitchState')
        val=logical(val);return
    end
catch
end
isLogical=false;
end
function varargout=var2str(varargin)
%Analogous to func2str, return the variable names as char arrays, as detected by inputname
% This returns an error for invalid inputs and if nargin~=max(1,nargout).
%
% You can use comma separated lists to create a cell array:
%   out=cell(1,2);
%   foo=1;bar=2;
%   [out{:}]=var2str(foo,bar);

% One-line alternative: function out=var2str(varargin),out=inputname(1);end
err_flag= nargin~=max(1,nargout) ;
if ~err_flag
    varargout=cell(nargin,1);
    for n=1:nargin
        try varargout{n}=inputname(n);catch,varargout{n}='';end
        if isempty(varargout{n}),err_flag=true;break,end
    end
end
if err_flag
    error('Invalid input and/or output.')
end
end