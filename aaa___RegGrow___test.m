function aaa___RegGrow___test(varargin)
%run a test suite

%these syntaxes test basic availability of functionalities and test cases
%that caused issues after writing the initial version of this function

if nargin==0,RunTestHeadless=false;else,RunTestHeadless=true;end

%% generate an example image
if exist('OCTAVE_VERSION', 'builtin')
    IM = get(0,'DefaultImageCData');IM = IM/max(IM(:));
    %Hash=ComputeNonCryptHash(RegGrow,256),return
    Hash='ED034137825C891268FA072122561C03D3404EF9963C2BDFCD60B9F9C93C3727';
    disp(['NOTE: different example image, ',...
        'performance can''t be compared with Matlab'])
else
    defimage = pow2(get(0,'DefaultImageCData'),47);
    IM = bitshift(defimage,-37);IM = fix(IM);
    IM = bitand(IM,31);IM = IM/max(IM(:));
    %Hash=ComputeNonCryptHash(RegGrow,256),return
    Hash='C2B3078720654204317BB886660926E63FADD194F67B03715D407593BD6EA236';
    
    % suppress v6.5 warning
    v=version;v(strfind(v,'.'):end)='';v=str2double(v);
    if v<=7
        warning off MATLAB:m_warning_end_without_block;clc
    end
end

%% run tests
clc
%attempt speedtest
try
    [sz,p]=speedtest;
    fprintf('speedtest finished\n\n')
    txt={sprintf('The waitbar for test 07 should close at %.0f%%\n',p)};
    txt{2}=sprintf('No waitbar should appear for test 08\n');
    if ~RunTestHeadless
        try%wrap in try-catch for compatibility with Octave 4.4.1
            uiwait(msgbox(txt))
        catch
        end
    end
catch
    sz=[50 50];
    fprintf('speedtest failed\n\n')
end
fails=0;
tic
try
    test01(Hash);
    fprintf('test 01: pass (%.1fs)\n',toc)
catch
    disp('test 01: fail')
    fails=fails+1;
end
drawnow;pause(0.01);tic
try
    test02(IM);
    fprintf('test 02: pass (%.1fs)\n',toc)
catch
    disp('test 02: fail')
    fails=fails+1;
end
drawnow;pause(0.01);tic
try
    test03(IM);
    fprintf('test 03: pass (%.1fs)\n',toc)
catch
    disp('test 03: fail')
    fails=fails+1;
end
drawnow;pause(0.01);tic
try
    test04(IM);
    fprintf('test 04: pass (%.1fs)\n',toc)
catch
    disp('test 04: fail')
    fails=fails+1;
end
drawnow;pause(0.01);tic
try
    test05(IM);
    fprintf('test 05: pass (%.1fs)\n',toc)
catch
    disp('test 05: fail')
    fails=fails+1;
end
drawnow;pause(0.01);tic
try
    test06(IM);
    fprintf('test 06: pass (%.1fs)\n',toc)
catch
    disp('test 06: fail')
    fails=fails+1;
end
drawnow;pause(0.01);tic
try
    test07(sz)
    fprintf('test 07: pass (%.1fs)\n',toc)
catch
    disp('test 07: fail')
    fails=fails+1;
end
drawnow;pause(0.01);tic
try
    test08(sz,p);
    fprintf('test 08: pass (%.1fs)\n',toc)
catch
    disp('test 08: fail')
    fails=fails+1;
end
drawnow;pause(0.01);tic
try
    test09(IM);
    fprintf('test 09: pass (%.1fs)\n',toc)
catch
    disp('test 09: fail')
    fails=fails+1;
end
drawnow;pause(0.01);tic
try
    test10(IM);
    fprintf('test 10: pass (%.1fs)\n',toc)
catch
    disp('test 10: fail')
    fails=fails+1;
end
if fails>0
    error('There were %d failed tests',fails)
else
    disp('test completed')
end
end
function [sz,p]=speedtest
%Determine the size needed to achieve a fixed runtime.
%This time should be enough to see at what approximate percentage the waitbar closes and should be
%long enough to notice that there is no waitbar opening.
%The second objective is to use the default inputs to determine the expected percentage for this
%kind of input.
%
%The required time probably doesn't scale linearly with the number of elements, but it should be
%close enough to get there within a handful of passes.
%(when tested across multiple releases and PCs, 3 or 4 passes seems to be the standard, with one
%outlier to 6 passes on the first run)

sz=[100 100];
target_time=2;%seconds
pass_count=0;
while true
    pass_count=pass_count+1;
    fprintf('running speedtest (pass %d) ',pass_count)
    old_elems=prod(sz);
    IM=1:prod(sz);
    IM=reshape(IM,sz);
    IM(:,2:2:end)=IM(end:-1:1,2:2:end);
    tic
    c=RegGrow(IM,'waitbar',true);
    testtime=toc;
    if abs(testtime-target_time)/target_time<0.2
        break
    end
    %adjust size to make this test last about 2 seconds
    elems=round(prod(sz)*target_time/testtime);
    sz=round(sqrt(elems))*[1 1];
    sz=2*ceil(sz/2);
    txt=sprintf('%.1E to %.1E',old_elems,prod(sz));
    txt=strrep(txt,'E+',' E');
    fprintf('(time=%.1f, elem count updated from %s)\n',testtime,txt)
end
p=100*sum(c(:))/prod(sz);
end
function test01(Hash)
try
    if ~strcmp(Hash,ComputeNonCryptHash(RegGrow))
        error('hash doesn''t match (did the default image change?)')
    end
    %pass
catch
    error('fail')
end
end
function test02(IM)
try
    opts=struct('seed',[6 4 2]);
    RegGrow(IM,opts);
    error('fail')
catch
    ME=lasterror; %#ok<LERR>
    if strcmp(ME.identifier,'HJW:RegGrow:incorrect_input_opt_seed')
        %pass
    else
        error('fail')
    end
end
end
function test03(IM)
try
    opts=struct('seed',[3;2]);
    RegGrow(IM,opts);
    %pass
catch
    error('fail')
end
end
function test04(IM)
try
    [out1,out2]=RegGrow(IM,0.3,[1 26],'overflow',1e2); %#ok<ASGLU>
    %pass
catch
    error('fail')
end
end
function test05(IM)
try
    [out1,OverflowFlag]=RegGrow(IM,'overflow',[]); %#ok<ASGLU>
    if OverflowFlag~=false
        error('fail')
    else
        %pass
    end
catch
    error('fail')
end
end
function test06(IM)
if exist('OCTAVE_VERSION', 'builtin')
    IM=repmat(IM,[1,1,5]);
    randfactor=20;
else
    IM=repmat(IM,[1,1,40]);
    randfactor=10;
end
randfactor=1+rand(size(IM))/randfactor;
IM=IM.*randfactor;
try
    [out1,OverflowFlag]=RegGrow(IM,'waitbar',false); %#ok<ASGLU>
    %pass
catch
    error('fail')
end
end
function test07(sz)
IM=1:prod(sz);
IM=reshape(IM,sz);
IM(:,2:2:end)=IM(end:-1:1,2:2:end);
try
    RegGrow(IM,'waitbar',true);
    %pass
catch
    error('fail')
end
end
function test08(sz,p)
IM=1:prod(sz);
IM=reshape(IM,sz);
IM(:,2:2:end)=IM(end:-1:1,2:2:end);
try
    %same test as 07, but now with the silent switch instead of waitbar
    result=RegGrow(IM,'silent',true);
    q=round(p*numel(result)/100);
    if abs( sum(result(:)) - q ) > 1
        error('fail')
    end
    %pass
catch
    error('fail')
end
end
function test09(IM)
try
    kernel=logical(zeros(3,3));%#ok<LOGL>
    w=warning('off','HJW:RegGrow:CentroidOnlyKernel');
    result=RegGrow(IM,'nHood',5,'kernel',kernel);%nhood is invalid, but is skipped
    warning(w);
    if sum(result(:))~=1
        error('fail')
    end
    %pass
catch
    error('fail')
end
end
function test10(IM)
try
    kernel=logical(zeros(3,3));%#ok<LOGL>
    kernel(end)=true;%the center voxel of the kernel doesn't actually matter
    result=RegGrow(IM,'kernel',kernel,'maxdiff',1+max(IM(:))-min(IM(:)));
    if sum(result(:))~=min(size(IM))
        error('fail')
    end
    %pass
catch
    error('fail')
end
end
function data=cast_to_uint16_vector(data)
%Linearize the input data and convert it to a uint16 vector
data=cast_to_uint16_vector__cell({data});
data([end-1 end])=[];%Remove the [1 1] that is added because of the wrapping in a cell
end
function data=cast_to_uint16_vector__cell(data)
sz=size(data);data=data(:);
for n=1:numel(data)
    if numel(data{n})==0
        c=double(class(data{n})');
        data{n}=uint16([0;c;size(data{n})']);
        continue
    end
    switch class(data{n})
        case {'double','single'}
            data{n}=cast_to_uint16_vector__floats(data{n});
        case 'logical'
            data{n}=cast_to_uint16_vector__logical(data{n});
        case {'uint8','uint16','uint32','uint64','int8','int16','int32','int64'}
            data{n}=cast_to_uint16_vector__integer(data{n});
        case 'char'
            data{n}=cast_to_uint16_vector__char(data{n});
        case 'cell'
            data{n}=cast_to_uint16_vector__cell(data{n});
        case 'struct'
            data{n}=cast_to_uint16_vector__struct(data{n});
        otherwise
            error('HJW:cast_to_uint16_vector:nosupport',...
                'Unsupported data type in nested variable')
    end
end
data=cell2mat(data);%Merge all cell contents
data=[data;uint16(mod(sz',2^16))];%Append the array size to make it influence the hash.
end
function data=cast_to_uint16_vector__floats(data)
sz=size(data);c=class(data);%the rest of the function treats singles as double

%Convert to a uint64, separate it into 4 words and everything merge into a vector
[bit64,bit4]=typecast_double_uint64(double(data));
bit4_round=mod(bit64,2^16);bit64=bit64-bit4_round;bit64=bit64/2^16;
bit3      =mod(bit64,2^16);bit64=bit64-bit3;      bit64=bit64/2^16;
bit2      =mod(bit64,2^16);bit64=bit64-bit2;      bit64=bit64/2^16;
bit1      =mod(bit64,2^16);
data=[bit1';bit2';bit3';bit4'];
data=uint16(data(:));

%Append the array size and type to make it influence the hash.
data=[data;uint16(c.');uint16(mod(sz',2^16))];
end
function data=cast_to_uint16_vector__logical(data)
sz=size(data);data=data(:);

if mod(numel(data),16) %pad to 16 bits with 0
    data(16*ceil(numel(data)/16))=0;
end
vector=uint16(2.^(15:-1:0))';
data=uint16(reshape(data,16,[]));
try
    data=data.*vector;
catch %no implicit expansion
    %times() is not defined for uint16 on ML6.5
    data=double(data).*repmat(double(vector),[1 size(data,2)]);
    data=uint16(data);
end
data=uint16(sum(data,1)).';

data=[data;uint16(mod(sz',2^16))];%Append the array size to make it influence the hash.
end
function data=cast_to_uint16_vector__integer(data)
%Large values (>2^52) will not have integer precision due to a conversion to double.
%This conversion is done, because of limited availability of operations on ML6.5.
sz=size(data);data=data(:);

c=class(data);
if c(1)~='u'
    %Shift int* values to the uint* data range
    data=double(data)-double(eval([c '(-inf)']));
else
    data=double(data);
end
switch c(end)
    case '8'
        %append a 0 for odd length and merge pairs
        if mod(numel(data),2),data(end+1)=0;end
        data=reshape(data,[],2);
        data=data(:,1)*255+data(:,2);
        data=uint16(data);
    case '6'
        data=uint16(data);
    case '2'
        %split to 2 words
        bit1=floor(data/2^16);
        bit2=mod(data,2^16);
        data=[bit1';bit2'];
        data=uint16(data(:));
    case '4'
        %split to 4 words, bit4 contains a rounding error for data >2^52
        bit64=data;
        bit4_round=mod(bit64,2^16);bit64=bit64-bit4_round;bit64=bit64/2^16;
        bit3      =mod(bit64,2^16);bit64=bit64-bit3;      bit64=bit64/2^16;
        bit2      =mod(bit64,2^16);bit64=bit64-bit2;      bit64=bit64/2^16;
        bit1      =mod(bit64,2^16);
        
        data=[bit1';bit2';bit3';bit4_round'];
        data=uint16(data(:));
end
%Append the array size and type to make it influence the hash.
data=[data;uint16(c.');uint16(mod(sz',2^16))];
end
function data=cast_to_uint16_vector__char(data)
sz=size(data);data=data(:);
data=uint16(data);%Normal chars are already a uint16 internally (in Matlab).
data=[data;uint16(mod(sz',2^16))];%Append the array size to make it influence the hash.
end
function data=cast_to_uint16_vector__struct(data)
sz=size(data);data=data(:);
fn=fieldnames(data);
output=cell(2,numel(fn));
for n=1:numel(fn)
    output{1,n}=fn{n};
    output{2,n}={data.(fn{n})};
end
data=cast_to_uint16_vector__cell(output);
data=[data;uint16(mod(sz',2^16))];%Append the array size to make it influence the hash.
end
function hash=ComputeNonCryptHash(data,varargin)
%Compute a non-cryptographic hash
%
% This function is intended to be fast, but without requiring a Java or mex implementation to do
% the actual hashing. It was *not* checked for any security flaws and is therefore probably
% vulnerable to most attacks.
% Non-cryptographic hashes should only be used as a checksum. Don't use this to do things like
% storing passwords.
%
%syntax:
%  hash=ComputeNonCryptHash(data)
%  hash=ComputeNonCryptHash(data,HashLength)
%
%data       The data to be hashed. Most common data types are allowed: uint*, int*, char, cell,
%           struct, double, or single (string is cast to char). The contents of the nested data
%           types (i.e. cell and struct) must also be one of the mentioned data types.
%HashLength The length of the hash (the number of bits). This value must be a multiple of 16. The
%           default is 256 bits. Depending on your input 64 bits might have some collisions, but
%           64 bits and higher should be safe.
%
%hash       The hash in an upper case hexadecimal char vector of size 1x(HashLength/4).
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
% Version: 1.0.1
% Date:    2020-07-06
% Author:  H.J. Wisselink
% Licence: CC by-nc-sa 4.0 ( https://creativecommons.org/licenses/by-nc-sa/4.0 )
% Email = 'h_j_wisselink*alumnus_utwente_nl';
% Real_email = regexprep(Email,{'*','_'},{'@','.'})

if nargin<1
    error('HJW:ComputeNonCryptHash:InputIncorrect','At least 1 input required.')
else
    try
        if isa(data,'string'),data=char(data);end
    catch
        error('HJW:ComputeNonCryptHash:ConvertToChar',...
            'The required conversion from string to char failed.')
    end
end
if nargin<2
    HashLength=256;
else
    HashLength=varargin{1};
    if numel(HashLength)~=1 || ~isnumeric(HashLength) || mod(HashLength,16)~=0 || HashLength<16
        error('HJW:ComputeNonCryptHash:InputIncorrect',...
            'Second input (hash length) must be a multiple of 16.')
    end
end

try
    %Convert the input to an uint16 array (Nx1).
    data=cast_to_uint16_vector(data);
catch
    ME=lasterror; %#ok<LERR>
    if strcmp(ME.identifier,'MATLAB:nomem')
        %rethrow memory error
        rethrow(ME)
    else
        error('HJW:ComputeNonCryptHash:UnwindFailed',...
            'The nested input contains an unsupported data type.')
    end
end

%Extend to a multiple of HashLength bits. Padding with zeros is generally not advised, and the
%performance penalty for this extension (compared to padding with zeros) should be negligible.
if mod(numel(data),HashLength/16)
    extra=uint16(1:HashLength/16).'; extra(1:mod(numel(data),HashLength/16))=[];
    data=[data;extra];
end

%Add perturbation to the data and convert to 16xN logical. Then further perturb the intermediate
%result by doing an XOR in blocks with itself (by reshaping and transposing).
data=ComputeNonCryptHash_shuffle_uint16(data);
data=ComputeNonCryptHash_uint16_to_logical(data);
data=xor(data,reshape(data,[],16).');

%Reshape to HashLength cols and collapse the key size down to the hash length by counting the
%number of true bits (even=1, odd=0).
data=mod(sum(reshape(data,HashLength,[]),2),2);
data=ComputeNonCryptHash_logical_to_uint16(data);

if nargin>2
    hash=data;%Return uint16 for the salting.
    return
end

%Perturb the hash, analogous to salting. This function computes the hash of the hash and applies a
%few operations to the data to increase the randomness of the end result.
data=ComputeNonCryptHash_add_salt(data);

%Convert the (HashLength/16)x1 uint16 to a hash string by encoding it as hexadecimal.
hash=ComputeNonCryptHash_dec2hex(data);hash=reshape(hash.',1,[]);
end
function data=ComputeNonCryptHash_add_salt(data)
%Apply a few transformations to the hash to increase the spread.
%If this function is not added, the hashes of -12345678 and 12345678 will be very similar.
%A true salt would be to append the salt to the data, so this is not actually a salt.
saltHashLength=16*numel(data);
salt=ComputeNonCryptHash(data,saltHashLength,[]);%avoid infinite recursion by using a third input
salt=ComputeNonCryptHash_shuffle_uint16_inv(salt);
data=mod(double(data).*double(salt),1+2^16);
data=uint16(data);
end
function hash=ComputeNonCryptHash_dec2hex(data)
%Look up the precomputed dec2hex for faster conversion.
persistent LUT
if isempty(LUT)
    LUT=upper(dec2hex(0:(-1+2^16),4));%even though the default should already upper case
end
data=double(data)+1;% plus() is not defined for uint16 on ML6.5
hash=LUT(data,:);
end
function data=ComputeNonCryptHash_logical_to_uint16(data)
if mod(numel(data),16) %pad to 16 bits with 0
    data(16*ceil(numel(data)/16))=0;
end
vector=uint16(2.^(15:-1:0))';
data=uint16(reshape(data,16,[]));
try
    data=data.*vector;
catch %no implicit expansion
    %times() is not defined for uint16 on ML6.5
    data=double(data).*repmat(double(vector),[1 size(data,2)]);
    data=uint16(data);
end
data=uint16(sum(data,1)).';
end
function data=ComputeNonCryptHash_shuffle_uint16(data)
%input should be uint16
base=65537;%base=(1+(2^16));
key=479001600;%key=1*2*3*4*5*6*7*8*9*10*11*12;
data = uint16(mod(double(data) * key , base));
end
function data=ComputeNonCryptHash_shuffle_uint16_inv(data)
base=65537;%base=(1+(2^16));
%key=1*2*3*4*5*6*7*8*9*10*11*12;
% %solution suggested by John D'Errico, https://www.mathworks.com/matlabcentral/answers/81859
% [G,C]=gcd(key,base);invKey=mod(C,base);
invKey=1919;
data=uint16(mod(double(data) * invKey,base));
end
function data=ComputeNonCryptHash_uint16_to_logical(data)
%uint16 Nx1 vector in, logical 16xN array out
persistent LUT
if isempty(LUT)
    LUT=dec2bin(0:(-1+2^16))=='1';
    LUT=LUT.';
end
data=double(data)+1;% plus() is not defined for uint16 on ML6.5
data=LUT(:,data);
end
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
function [bit64,last16bits]=typecast_double_uint64(FP)
%Turn a double into a uint64 with the same binary representation.
%This is similar to typecast(FP,'uint64'); the difference being that this function ignores
%endianness, supports array inputs, and is slower.
%
%Because of missing support for some operations in ML6.5, this function returns the uint64 as a
%double. Because this may cause rounding errors, the last 16 bits are returned separately.

[M,E]=log2(FP);
signBit =-floor(sign(FP)/2-0.5);
exponent=E+1022;
mantissa=abs(M)*2-1;

%no plus() for integer types in ML6.5, so we need to use double, instead of uint64
bit64=zeros(size(FP));
bit64=bit64+(signBit*2^63);
bit64=bit64+(exponent*2^52);
bit64=bit64+(mantissa*2^52);
last16bits=mod(mantissa*2^52,2^16);

%correct the 0 and hard-code the special cases
L=isinf(FP);
bit64(FP==0)=0;
bit64(isnan(FP))=18444492273895866368;
bit64(L & FP>0)=9218868437227405312;%positive inf
bit64(L & FP<0)=18442240474082181120;%negative inf
last16bits(FP==0)=0;
last16bits(isnan(FP))=0;
last16bits(L)=0;
end