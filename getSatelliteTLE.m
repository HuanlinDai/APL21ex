function [TLE] = getSatelliteTLE(ID,inEpochDatenum)
%% *purpose*
% return the TLE for Satellite based on epoch
%% *inputs*
%  ID - Spacecraft ID
%       73027 = Skylab
%  inEpochDatenum - TLEs can be selected based on epoch 
%% *outputs*
%  TLE - the two line element set corresponding to the satellite at that
%        epoch
%% *history*
%  When       Who    What
%  ---------- ------ --------------------------------------------------
%  2019/07/17 mnoah  original code
%  2020/01/19 mnoah  placeholder

if (ID == 25544)
    TLE = { ...
        'ISS (ZARYA)'; ...
        '1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927'; ...
        '2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537'};
elseif (ID == 255440)
    TLE = { ...
        'ISS (ZARYA)' ...
        '1 25544U 98067A   04236.56031392  .00020137  00000-0  16538-3 0  9993' ...
        '2 25544  51.6335 344.7760 0007976 126.2523 325.9359 15.70406856328906'};
elseif (ID == 255442)
    TLE = { ...
        'ISS (ZARYA)' ...
        '1 25544U 98067A   21188.62787601  .00000229  00000-0  12373-4 0  9990' ...
        '2 25544  51.6436 241.6415 0001986 141.2144 354.6717 15.48768433291768'};
elseif (ID == 07376)
    TLE = { ...
        'MOLNIYA 2-10' ...         
        '1 07376U 74056A   21200.78494011 -.00000835  00000-0  00000-0 0  9991' ...
        '2 07376  63.7488 128.6248 6670497 289.4096  13.1550  2.01178017343663'};
elseif (ID == 13070)
    TLE = { ...
        'MOLNIYA 1-53' ...
        '1 13070U 82015A   21180.17498817 -.00001367  00000-0  52409-2 0  9994' ...
        '2 13070  62.4656 261.6018 7274093 280.5474  11.8756  2.00521228288163'};
elseif (ID == 07392)
    TLE = { ...
        'MOLNIYA 1-S' ...
        '1 07392U 74060A   21179.91528292 -.00000134  00000-0  00000+0 0  9994' ...
        '2 07392   5.9914 288.2114 0010874 194.2146 201.6455  1.00345146115963'};
elseif (ID == 08195)
    TLE = { ...
        'MOLNIYA 2-14' ...
        '1 08195U 75081A   21177.87857582  .00001190  00000-0  00000-0 0  9993' ...
        '2 08195  63.3181 228.7209 7041687 294.0542  10.0704  2.02060154335706'};
elseif (ID == 08601)
    TLE = { ...
        'MOLNIYA 1-32' ...
        '1 08601U 76006A   21186.57570768  .00001350  00000-0  14533-2 0  9994' ...
        '2 08601  63.1753 349.0945 7395835 282.2691  10.5046  1.99973651 45731'};
elseif(ID == 09880)
    TLE = { ...
        'MOLNIYA 1-36' ...
        '1 09880U 77021A   21186.31736846 -.00001812  00000-0  26004-3 0  9996' ...
        '2 09880  63.1385 318.8663 7384958 288.2227   9.3553  2.00769617222406'};
elseif(ID == 10455)
    TLE = { ...
        'MOLNIYA 3-8' ...
        '1 10455U 77105A   21186.62369222 -.00001145  00000-0  29403-3 0  9998' ...
        '2 10455  62.8200  52.3393 7380177 279.9913  11.1468  2.00800761320104'};
elseif(ID == 13875)
    TLE = { ...
        'MOLNIYA 3-20' ...
        '1 13875U 83015A   21186.22899311  .00001261  00000-0  58007-2 0  9996' ...
        '2 13875  63.2427  46.5921 7335821 256.0644  19.5571  2.00586652216745'};
elseif(ID == 13890)
    TLE = { ...
        'MOLNIYA 1-56'            
        '1 13890U 83019A   21186.01188897 -.00001019  00000-0  15357-2 0  9997' ...
        '2 13890  63.1315  72.3784 7341171 261.6580  17.1236  1.99953138279652'};
elseif(ID == 48274)
    TLE = { ...
        'TIANHE-1' ...
        '1 48274U 21035A   21188.18655527 -.00006823  00000-0 -70895-4 0  9999' ...
        '2 48274  41.4697 161.0904 0006875  77.6216 357.0670 15.62384497 10822'};
elseif (ID == 255443)
    TLE = { ...
        'ISS (ZARYA)' ...      
        '1 25544U 98067A   21197.40752575  .00002080  00000-0  46171-4 0  9998'...
        '2 25544  51.6421 198.2252 0001949 161.2303 359.1523 15.48808816293123'};
elseif (ID == 48864)
    TLE = { ...
        'TIANQI-14'...               
        '1 48864U 21055E   21197.16701608 -.00000298  00000-0 -26550-4 0  9995' ...
        '2 48864  35.0046 172.8401 0007083 165.2555 194.8355 14.91395958  3866'};
elseif (ID == 37849)
    TLE = { ...
        'SUOMI NPP' ...              
        '1 37849U 11061A   21197.36024387 -.00000021  00000-0  11037-4 0  9993'...
        '2 37849  98.7387 135.6491 0000875 105.0440  54.3027 14.19536906503506'};
elseif (ID == 19548)
    TLE = { ...
        'TDRS 3' ...                 
        '1 19548U 88091B   21201.78210319 -.00000307  00000-0  00000-0 0  9991' ...
        '2 19548  13.9418 353.9933 0033981 327.7475 209.7980  1.00271450107417'};
elseif (ID == 3002)
    TLE = { ...
        'M00N 1' ...
        '1 00000U 01001A   21202.79722222 -.00000298  00000-0  00000-4 0  9995' ...
        '2 00000  18.2800 318.1500 0055400 165.2555 135.2700 00.03660000999999'};
else
    error('ID not found or coded yet; modify code to include more satellites');
end

end




