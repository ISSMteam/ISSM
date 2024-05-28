%
%  convert the units for dakota responses.
%
%  [dresp]=dakota_resp_uconv(dresp)
%
%  where the required input is:
%    dresp         (structure array, responses)
%
%  the required output is:
%    dresp         (structure array, responses)
%
%  this function reads through a dakota response structure, and
%  for those quantities whose descriptors are recognized, converts
%  the units of all the applicable fields.  a "unit" field is added
%  to the response structure.
%
%  this data would typically be read by dakota_out_parse and be used
%  for plotting and other post-processing within matlab or excel.
%
%  "Copyright 2010, by the California Institute of Technology.
%  ALL RIGHTS RESERVED. United States Government Sponsorship
%  acknowledged. Any commercial use must be negotiated with
%  the Office of Technology Transfer at the California Institute
%  of Technology.  (NTR 47078)
%
%  This software may be subject to U.S. export control laws.
%  By accepting this  software, the user agrees to comply with
%  all applicable U.S. export laws and regulations. User has the
%  responsibility to obtain export licenses, or other export
%  authority as may be required before exporting such information
%  to foreign countries or providing access to foreign persons."
%
function [dresp]=dakota_resp_uconv(dresp)

if ~nargin
    help dakota_resp_uconv
    return
end

if ~isstruct(dresp)
    error('''%s'' is not a structure array.',inputname(1));
end
if ~isfield(dresp,'descriptor')
    error('''%s'' does not have a descriptor field.',inputname(1));
end

%%  define the conversion factors

sec_per_yr=365.2425*24*60*60;    %  mean gregorian year
m_per_km=1000;
kg_per_gton=10^12;

%%  loop through the response array

for i=1:numel(dresp)
    dresp(i).unit='';
    if     ~isempty(strfind(dresp(i).descriptor,'vel'))    %  in m/sec
        dresp(i)=drespi_conv(dresp(i),sec_per_yr,'m/yr');
    elseif ~isempty(strfind(dresp(i).descriptor,'misfit'))    %  in m^2*(m/sec)^2
        dresp(i)=drespi_conv(dresp(i),1/m_per_km^2*sec_per_yr^2,'km^2*(m/yr)^2');
    elseif ~isempty(strfind(dresp(i).descriptor,'mass_flux'))    %  in kg/sec
        dresp(i)=drespi_conv(dresp(i),1/kg_per_gton*sec_per_yr,'Gton/yr');
    else
        disp(['Skipping response ''' dresp(i).descriptor '''.']);
    end
end

end

%%  function to convert the units of a dakota response

function [dresp]=drespi_conv(dresp,unew_per_uold,ulab)

disp(['Converting response ''' dresp.descriptor ''' to ' ulab '.']);

%  loop over the fields, converting only the appropriate ones

fnames=fieldnames(dresp);
for i=1:length(fnames)
    switch fnames{i}
        case {'sample',...
              'mean',...
              'stddev',...
              'meanci',...
              'stddevci',...
              'min',...
              'quart1',...
              'median',...
              'quart3',...
              'max'}    %  appropriate to convert
            dresp.(fnames{i})=dresp.(fnames{i})*unew_per_uold;
        case {'cdf'}    %  only responses, not probs or reliabilities
            dresp.cdf(:,1)=dresp.cdf(:,1)*unew_per_uold;
        case {'descriptor',...
              'coefvar',...
              'var',...
              'impfac'}    %  unitless (or non-numeric)
            continue;
        case {'best',...
              'vum',...
              'unit'}    %  other
            continue;
        otherwise
            disp(['Unrecognized field ''' fnames{i} '''.']);
    end
end

if isfield(dresp,'unit')
    dresp.unit=ulab;
end

end
