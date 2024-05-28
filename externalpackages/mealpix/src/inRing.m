function  pixList = inRing(nSide, nRing, varargin)
% INRING return pixels in ring or ring section
%
%  pixList = inRing(nSide, nRing, phi0, dPhi, 'Param1', Value1, ...);
%
%  nSide       HEALPix resolution parameter
%  nRing       Ring number (1 <= nRing <= 4*nSide-1)
%  phi0, dPhi  (optional) Ring section longitudes [phi-dPhi,phi+dPhi]
%
%  Param       Value
%  'nest'      return pixels with nested indexing (true | {false})
%  'strict'    return only pixels whose center is in longitude range
%              (true | {false})
%
%  pixList     size(nRing) cell array of MEALPix pixels in ring defined by
%              nRing, phi0, dPhi
%
%  nRing may be a numeric array, in which case phi0, dPhi and nest may each
%  be either scalar or a size(nRing) array. If strict is false (default)
%  then all pixels intersected by longitude range are included in pixList.
%
% Authors: Lee Samuel Finn, Matthew Kinsey
% Copyright 2010-2011 Lee Samuel Finn

% $Id: inRing.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

% Copyright 2010-2011 by Lee Samuel Finn. 
% This program file is part of MEALPix. It is licensed under the Apache
% License, Version 2.0 (the  "License"); you may not use MEALPix except in
% compliance with the License. You may obtain a copy of the License at 
% <http://www.apache.org/licenses/LICENSE-2.0> 
%
% Unless required by applicable law or agreed to in writing, MEALPix
% software distributed under the License is distributed on an "AS IS"
% BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
% implied. See the License for the specific language governing permissions
% and limitations under the License.


%% Parse and validate arguments
p = inputParser;
p.addRequired('nSide',@qnSide);
p.addRequired('nRing',@(n)(qNRing(nSide,n)));
p.addOptional('phi0',0,@isnumeric);
p.addOptional('dPhi',2*pi,@isnumeric);
p.addParamValue('nest',false,@islogical);
p.addParamValue('strict',false,@islogical);

p.parse(nSide,nRing,varargin{:});

phi0 = p.Results.phi0;
dPhi = p.Results.dPhi;
nest = p.Results.nest;
strict = p.Results.strict;

% Validate dimensions
szRing = size(nRing);

if 1 == prod(szRing)
  if all(1 ~= [numel(phi0), numel(dPhi), numel(nest), numel(strict)])
    msgid = [mfilename ':dims'];
    error(msgid,msgid);
  end
else
  if 1 < numel(phi0)
    assert(all(size(phi0) == szRing),[mfilename ':phi0Dims']);
  else
    phi0 = repmat(phi0,szRing);
  end
  if 1 < numel(dPhi)
    assert(all(size(dPhi) == szRing),[mfilename ':dPhiDims']);
  else
    dPhi = repmat(dPhi,szRing);
  end
  if 1 < numel(nest)
    assert(all(size(nest) == szRing),[mfilename ':nestDims']);
  else
    nest = repmat(nest,szRing);
  end
  if 1 < numel(strict)
    assert(all(size(strict) == szRing),[mfilename ':strictDims']);
  else
    strict = repmat(strict,szRing);
  end
end

%% Initializations
npix = nSide2nPix(nSide);
ncap  = 2*nSide*(nSide-1);
pixList = cell(size(nRing));

%% Put phi0, dPhi in canonical form
dPhi = abs(dPhi);
phi_low = mod(phi0 - dPhi,2*pi);
phi_hi  = mod(phi0 + dPhi,2*pi);
mLong = dPhi >= pi;
if 1 == numel(mLong)
  mLong = repmat(mLong,size(nRing));
end

%% Identify ring number
nr = zeros(size(nRing));       % ring number
iPix1 = zeros(size(nRing));    % lowest pixel number in the ring
iPix2 = zeros(size(nRing));    % highest pixel number in the ring
kShift = zeros(size(nRing));

%%% equatorial region rings
mask = nSide <= nRing & nRing <= 3*nSide;
mir = nRing(mask) - nSide + 1;
iPix1(mask) = ncap  + 4*nSide*(mir-1);
iPix2(mask) = iPix1(mask) + 4*nSide - 1;
kShift(mask)  = mod(mir,2);
nr(mask) = 4*nSide;

%%% north pole rings
mask = nRing < nSide;
mir = nRing(mask);
iPix1(mask) = 2*mir.*(mir-1);
iPix2(mask) = iPix1(mask) + 4*mir - 1;
kShift(mask) = 1;
nr(mask) = 4*mir;

%%% south pole rings
mask = nRing > 3*nSide;
mir = 4*nSide - nRing(mask);
iPix1(mask) = npix - 2*mir.*(mir+1);
iPix2(mask) = iPix1(mask) + 4*mir - 1;
kShift(mask) = 1;
nr(mask) = 4*mir;

%% Construct pixel list: Full ring case
if any(mLong(:))
  %%% Save all ring pixels: ring indexing
  temp = arrayfun(...
    @(i1,i2)(i1:i2),...
    iPix1(mLong),...
    iPix2(mLong),...
    'UniformOutput',false);
  [pixList{mLong}] = deal(temp{:});
end

%% Construct pixel list: Partial ring case
if any(~mLong(:))
  shift = kShift/2;
  
  nip = sum(~mLong(:));
  ip_low = zeros(nip,1);
  ip_hi = zeros(nip,1);
  nList = false(nip,1);
  
  if any(~strict(:))
    % ~strict = conservative : include every intersected pixel,
    % even if pixel CENTER is not in the range [phi_low, phi_hi]
    mask = ~strict & ~mLong;
    mnr = nr(mask);
    mip_low = round(mnr.*phi_low(mask)/(2*pi) - shift(mask));
    mip_hi  = round(mnr.*phi_hi(mask)/(2*pi) - shift(mask));
    ip_low(~strict(~mLong)) = mod(mip_low, mnr);% in {0,nr-1}
    ip_hi(~strict(~mLong))  = mod(mip_hi, mnr); % in {0,nr-1}
  end
  
  if any(strict)
    mask = strict & ~mLong;
    mnr = nr(mask);
    % strict : include only pixels whose CENTER is in [phi_low, phi_hi]
    mip_low = ceil(mnr.*phi_low(mask)/(2*pi) - shift(mask));
    mip_hi  = floor(mnr.*phi_hi(mask)/(2*pi) - shift(mask));
    df = mod(mip_low - mip_hi, mnr);
    df(df<0) = df(df<0) + mnr(df<0);
    
    % Mark intervals so small (and away from pixel center) that no pixel
    % is included in it
    n0 = (df == 1) & (dPhi(mask).*mnr < pi);
    nList(n0) = true;
    
    m0 = mip_low >= mnr;
    mip_low(m0) = mip_low(m0) - mnr(m0);
    m0 = mip_hi < 0;
    mip_hi(m0) = mip_hi(m0) + mnr(m0);
    
    ip_low(strict(~mLong)) = mip_low;
    ip_hi(strict(~mLong)) = mip_hi;
  end;
  
  ip_low = ip_low + iPix1(~mLong(:));
  ip_hi  = ip_hi  + iPix1(~mLong(:));
  
  % identify rings with "wrap-around" longitudinal slices
  mTop = ip_low > ip_hi;
  ip2 = iPix2(~mLong(:));
  ip1 = iPix1(~mLong(:));
  
  topList = arrayfun(...
    @(pLow,p1,p2,pHi)([pLow:p2, p1:pHi]),...
    ip_low(mTop),ip1(mTop),ip2(mTop),ip_hi(mTop), ...
    'UniformOutput',false);
  
  % rings with simply-connected longitudinal slices
  nTopList = arrayfun(...
    @(l,h)(l:h),...
    ip_low(~mTop),ip_hi(~mTop),...
    'UniformOutput',false);
  
  mList = cell(sum(~mLong(:)),1);
  [mList{mTop}] = deal(topList{:});
  [mList{~mTop}] = deal(nTopList{:});
  [mList{nList}] = deal([]);
  
  [pixList{~mLong}] = deal(mList{:});
  
end

%% Convert ring indexed to nest index as required
if any(nest)
  [nestList{1:sum(nest(:))}] = deal(pixList{nest});
  nestList = cellfun(...
    @(n)(ring2nest(nSide,n)),...
    nestList, 'UniformOutput',false);
  [pixList{nest}] = deal(nestList{:});
end

%% Convert to MEALPix indexing
pixList = cellfun(@(x)(x+1),pixList,'UniformOutput',false);

return

function tf = qNRing(nSide,nRing)
tf = all(1 <= nRing(:)) && all(nRing(:) < 4*nSide);
return
