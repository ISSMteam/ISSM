function [c] = geoCode(address, service, key)
%GEOCODE look up the latitude and longitude of a an address
%
%   COORDS = GEOCODE( ADDRESS ) returns the geocoded latitude and longitude 
%   of the input address. 
%
%   COORDS = GEOCODE( ADDRESS, SERVICE) performs the look up using the
%   specified SERVICE. Valid services are
%       google  - Google Maps  (default service)
%       osm     - OpenStreetMap
%       yahoo   - Yahoo! Place Finder 
%
%   COORDS = GEOCODE( ..., SERVICE, APIKEY) allows the specifcation of an AppId
%   API key if needed.

% Copyright(c) 2012, Stuart P. Layton <stuart.layton@gmail.com>
% https://stuartlayton.com
%
% Revision History
%   2012/08/20 - Initial Release
%   2012/08/20 - Simplified XML parsing code


% Validate the input arguments

% Check to see if address is a valid string
if isempty(address) || ~ischar(address) || ~isvector(address)
    error('Invalid address provided, must be a string');
end

% if no service is specified or an empty service is specified use google
if nargin<2 || isempty(service) 
    service = 'google';
end

% if no key is specified then set it to empty, also check to see if char array
if nargin<3 
    key = [];
end

%replace white spaces in the address with '+'
address = regexprep(address, ' ', '+');

% Switch on the specified service, construct the Query URL, and specify the
% function that will be used to parse the resulting XML 
switch lower(service)
    case('google')
        SERVER_URL = 'https://maps.google.com';
        queryUrl = sprintf('%s/maps/geo?output=xml&q=%s',SERVER_URL, address);
        parseFcn = @parseGoogleMapsXML;

    case('yahoo')
      
        SERVER_URL = 'https://where.yahooapis.com/geocode';
        queryUrl = sprintf('%s?location=%s',SERVER_URL, address);
       
% The Yahoo docs say that an AppID is required although
% it appears that responses are given without a valid appid 
% If an AppId is provided include it in the URL
        if ~isempty(key)
            queryUrl = sprintf('%s&appid=%s', queryUrl, key);
        end
        
        parseFcn = @parseYahooLocalXML;
        
    
    case {'osm', 'openstreetmaps', 'open street maps'}
        
        SERVER_URL = 'https://nominatim.openstreetmap.org/search';
        queryUrl = sprintf('%s?format=xml&q=%s', SERVER_URL, address);
        parseFcn = @parseOpenStreetMapXML;

    otherwise
        error('Invalid geocoding service specified:%s', service);
end

    try
        docNode = xmlread(queryUrl);
    catch  %#ok<CTCH>
        error('Error, could not reach %s, is it a valid URL?', SERVER_URL);
    end
   
    c = parseFcn(docNode);
end

% Function to parse the XML response from Google Maps
function [c] = parseGoogleMapsXML(docNode)
    
    %check the response code to see if we got a valid response
    codeEl = docNode.getElementsByTagName('code');
    errCode = str2double( char( codeEl.item(0).getTextContent ) );
    % code 200 is associated with a valid geocode xml file
    if errCode~=200
        fprintf('No data received from server! Received code:%d\n', errCode)
        c = nan(2,1);
        return;
    end

    %get the 'coordinates' element from the document
    cordEl = docNode.getElementsByTagName('coordinates');
    
    %make sure the xml actually included a coordinates tag
    if cordEl.length<1
        c = nan(2,1);
        warning('No coordinates returned for the specified address');
        return;
    end
        
    % get the coordinates from the first node, convert them to numbers
    coords = cellfun(@str2double, regexp( char( cordEl.item(0).getTextContent), ',', 'split'));
    
    c = coords([2, 1]); % return the latitude and longitude
end

% Function to parse the XML response from Yahoo Local
function [c] = parseYahooLocalXML(docNode)
    
    %check the response code to see if we got a valid response
    codeEl = docNode.getElementsByTagName('Error');
    errCode = str2double( char( codeEl.item(0).getTextContent ) );
    
    % code 0 is associated with a valid geocode xml file
    if errCode~=0
        fprintf('No data received from server! Received code:%d\n', errCode)
        c = nan(2,1);
        return;
    end
    
    %check to see if a location was actually found
    foundEl = docNode.getElementsByTagName('Found');
    found = str2double( char( foundEl.item(0).getTextContent) );
    % 
    if found<1
        disp('A location with that address was not found!');
        c = nan(2,1);
        return;
    end
    
    
    latEl = docNode.getElementsByTagName('latitude');
    lonEl = docNode.getElementsByTagName('longitude');
    
    %make sure the xml actually included latitude and longitude tags
    if latEl.length==0 || lonEl.length==0
        c = nan(2,1);
        disp('No coordinates were found for that address');   
        return;
    end
    
    c(1) = str2double( char( latEl.item(0).getTextContent) );
    c(2) = str2double( char( lonEl.item(0).getTextContent) );
    
end

% Function to parse the XML response from OpenStreetMap
function [c] = parseOpenStreetMapXML(docNode)
    
    serverResponse = docNode.getElementsByTagName('searchresults').item(0);
    placeTag = serverResponse.getElementsByTagName('place').item(0);
    
    if isempty(placeTag)
        disp('OpenStreeMap returned no data for that address');%#ok
        c = nan(2,1);
        return;
    end
    
    c(1) = str2double( char( placeTag.getAttribute('lat') ) );
    c(2) = str2double( char( placeTag.getAttribute('lon') ) );
    
end


function elementText = GetElementText(resultNode,elementName)
% GETELEMENTTEXT given a result node and an element name
% returns the text within that node as a Matlab CHAR array

elementText = ...
    char( resultNode.getElementsByTagName(elementName).item(0).getTextContent );

end
