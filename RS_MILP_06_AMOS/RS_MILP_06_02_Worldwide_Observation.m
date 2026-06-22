% ============================================================
% 54 Global Major City Ground Targets
% ============================================================

cityData = [
    {"Seoul",           37.57, 126.98}
    {"Busan",           35.18, 129.08}
    {"Incheon",         37.46, 126.71}
    {"Tokyo",           35.68, 139.69}
    {"Osaka",           34.69, 135.50}
    {"Nagoya",          35.18, 136.91}
    {"Beijing",         39.90, 116.41}
    {"Shanghai",        31.23, 121.47}
    {"Shenzhen",        22.54, 114.06}
    {"Hong Kong",       22.32, 114.17}
    {"Taipei",          25.03, 121.57}
    {"Singapore",        1.35, 103.82}

    {"Bangkok",         13.76, 100.50}
    {"Ho Chi Minh City",10.82, 106.63}
    {"Jakarta",         -6.21, 106.85}
    {"Manila",          14.60, 120.98}
    {"Kuala Lumpur",     3.14, 101.69}
    {"Mumbai",          19.08, 72.88}
    {"Delhi",           28.61, 77.21}
    {"Bengaluru",       12.97, 77.59}

    {"London",          51.51, -0.13}
    {"Paris",           48.86, 2.35}
    {"Amsterdam",       52.37, 4.90}
    {"Brussels",        50.85, 4.35}
    {"Frankfurt",       50.11, 8.68}
    {"Berlin",          52.52, 13.40}
    {"Munich",          48.14, 11.58}
    {"Madrid",          40.42, -3.70}
    {"Barcelona",       41.39, 2.17}
    {"Rome",            41.90, 12.50}
    {"Milan",           45.46, 9.19}
    {"Warsaw",          52.23, 21.01}

    {"New York",        40.71, -74.01}
    {"Washington DC",   38.90, -77.04}
    {"Boston",          42.36, -71.06}
    {"Chicago",         41.88, -87.63}
    {"Dallas",          32.78, -96.80}
    {"Houston",         29.76, -95.37}
    {"Los Angeles",     34.05,-118.24}
    {"San Francisco",   37.77,-122.42}
    {"Seattle",         47.61,-122.33}
    {"Toronto",         43.65, -79.38}
    {"Montreal",        45.50, -73.57}
    {"Mexico City",     19.43, -99.13}
    {"Austin",          30.26, -97.7431}

    {"Dubai",           25.20, 55.27}
    {"Riyadh",          24.71, 46.67}
    {"Tel Aviv",        32.09, 34.78}
    {"Cairo",           30.04, 31.24}
    {"Johannesburg",   -26.20, 28.04}
    {"Lagos",            6.52, 3.38}

    {"Sao Paulo",      -23.55,-46.63}
    {"Rio de Janeiro", -22.91,-43.17}
    {"Buenos Aires",   -34.60,-58.38}
    {"Santiago",       -33.45,-70.67}
];


cityNames = string(cityData(:,1));
lats      = cell2mat(cityData(:,2));
lons      = cell2mat(cityData(:,3));

% ============================================================
% Ground Station Generation
% ============================================================

N = numel(lats);

fprintf("Total %4.0f ground locations are generated\n", N)

gs = groundStation(sc, lats, lons, ...
    "Name", "Ground_Point_" + string(1:N)');

set(gs,'ShowLabel',0)
set(gs,'MarkerSize',3)
set(gs,'MarkerColor','m')


figure;
geoplot(lats, lons, ...
    'Marker','o', ...
    'LineStyle','none', ...
    'MarkerSize',5, ...
    'MarkerFaceColor','r')

geobasemap streets

geolimits([-40 60], [-130 150])

title('Global Ground Observation Points', ...
    'FontSize',12, ...
    'FontWeight','bold')

for k = 1:N
    text(lats(k), lons(k)+1, cityNames(k),'FontSize',8)
end