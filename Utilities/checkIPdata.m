function IPdata = checkIPdata
% checks IP data using ip-api.com and returns it as a Matlab structure
% see http://ip-api.com/docs/api:returned_values for more info

try
    x = (urlread('http://ip-api.com/csv'));
    x = regexp(x,',','split');
    IPdata.status = x{1};
    IPdata.country = x{2};
    IPdata.countryCode = x{3};
    IPdata.region = x{4};
    IPdata.regionName = x{5};
    IPdata.city = x{6};
    IPdata.zip = x{7};
    IPdata.lat = x{8};
    IPdata.lon = x{9};
    IPdata.timezone = x{10};
    IPdata.isp = x{11};
    IPdata.org = x{12};
    % IPdata.as = x{13};
    % IPdata.reverse = x{14};
    % IPdata.mobile = x{15};
    % IPdata.proxy = x{16};
    % IPdata.query = x{17};
    
catch
    IPdata = 'IP data unavailable';
end
