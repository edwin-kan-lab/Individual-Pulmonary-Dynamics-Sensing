function datetimeObj = labviewDateToDateTime(labviewStr, year)

month = labviewStr(1:2);
date = labviewStr(3:4);
hour = labviewStr(6:7);
min = labviewStr(8:9);
sec = labviewStr(10:11);
datetimeObj = datetime([year '-' month '-' date ' ' hour ':' min ':' sec], 'InputFormat', 'yyyy-MM-dd HH:mm:ss');


end