library(phylin);
library(lubridate);
library(raster);
library(rgdal);
library(tmap);
library(maptools);
library(spatstat);
library(gstat); 
library(sp);    
library(geosphere);
library(dplyr);

working_directory = "F:/Github/Analisis-Penertiban-Tanah-Abang/";
working_directory2 = "C:/Users/Fox/Documents/GitHub/Analisis-Penertiban-Tanah-Abang/";

setwd(working_directory2);

pasar_tanahabang = c(106.814535,-6.187611);

alerts = dbGetQuery(conn,"select * from sor.waze_alerts where date_id>=20171001 and city in ('Jakarta Pusat')");
jams = dbGetQuery(conn,"select distinct street,roadtype,city,line_x,line_y,delay,pubmillis from (select distinct street,roadtype,city,line_x,line_y,delay,pubmillis from sor.waze_jams where date_id>=20171028 and city in ('Jakarta Pusat')) t1");

alerts_29 = dbGetQuery(conn,"select * from sor.waze_alerts where date_id>20171028 and city in ('Jakarta Pusat')");
jams_29 = dbGetQuery(conn,"select distinct street,roadtype,city,line_x,line_y,delay,pubmillis from (select distinct street,roadtype,city,line_x,line_y,delay,pubmillis from sor.waze_jams where date_id>20171028 and city in ('Jakarta Pusat')) t1");

alerts$line_x = as.numeric(alerts$waze_alerts.location_x);
alerts$line_y = as.numeric(alerts$waze_alerts.location_y);
alerts$pubmillis = as.numeric(alerts$waze_alerts.pubmillis);

jams$line_x = as.numeric(jams$line_x);
jams$line_y = as.numeric(jams$line_y);
jams$pubmillis = as.numeric(jams$pubmillis);

jams = jams %>%
rowwise() %>%
mutate(distance=distGeo(c(line_x,line_y),pasar_tanahabang));

alerts = alerts %>%
rowwise() %>%
mutate(distance=distGeo(c(line_x,line_y),pasar_tanahabang));

#####################################################
## CALCULATE DISTANCE FROM EPICENTRUM
#####################################################

library(doParallel);
library(foreach);

calc_distance_parallel = function(alerts,loc){

cores = detectCores(logical = FALSE)-1;
cl = makeCluster(cores);
registerDoParallel(cl, cores=cores);
 
chunk.size = ceiling(nrow(alerts)/cores);
total.chunk.size = cores * chunk.size;
diff.chunk = total.chunk.size - nrow(alerts);

added_data = 0;

if(diff.chunk>0){
	for(i in 1:diff.chunk){
		alerts = rbind(alerts,alerts[1,]);
		added_data=added_data+1;
	}
}

start_time = Sys.time();

 res2.p = foreach(i=1:cores, .combine='rbind', .packages='geosphere') %dopar%
 { 
    res = matrix(0, nrow=chunk.size, ncol=1)
    for(x in ((i-1)*chunk.size+1):(i*chunk.size)) {
        res[x - (i-1)*chunk.size,] = distGeo(c(alerts$line_x[x],alerts$line_y[x]),loc);
    }
    res;
 }

end_time = Sys.time();
end_time - start_time;

if(diff.chunk>0){
	for(i in 1:diff.chunk){
		res2.p = res2.p[-nrow(res2.p),];
	}
}

alerts$jarak = res2.p[,1];

stopImplicitCluster();
stopCluster(cl);

return(alerts);

}

alerts_dist = calc_distance_parallel(alerts,pasar_tanahabang);
jams_dist = calc_distance_parallel(jams,pasar_tanahabang);

#####################################################
## FILTERING RADIUS 1 KM
#####################################################

jams_distance = calc_distance_parallel(jams,loc);
alerts_distance = calc_distance_parallel(alerts,loc);

jams_1km = jams_dist %>%
filter(jarak<=1000) %>%
as.data.frame();

alerts_1km = alerts_dist %>%
filter(jarak<=1000) %>%
as.data.frame();

library(lubridate);

jams_1km$Waktu = as_datetime(jams_1km$pubmillis/1000,tz=Sys.timezone());
jams_1km$weekday = wday(jams_1km$Waktu,label=TRUE);
jams_1km$hour = hour(jams_1km$Waktu);
jams_1km$yday = yday(jams_1km$Waktu);

alerts_1km$Waktu = as_datetime(alerts_1km$pubmillis/1000,tz=Sys.timezone());
alerts_1km$weekday = wday(alerts_1km$Waktu,label=TRUE);
alerts_1km$hour = hour(alerts_1km$Waktu);
alerts_1km$yday = yday(alerts_1km$Waktu);
alerts_1km$week = week(alerts_1km$Waktu);

alerts_1km$waktu = as.character(alerts_1km$Waktu);

#####################################################
## CALCULATE KERNEL DENSITY
#####################################################

load("alerts_jams_20171001_tanahabang.RData");

subtype = unique(alerts_1km$waze_alerts.subtype);
yday = unique(alerts_1km$yday)
week = unique(alerts_1km$week)

library(sp)
e = as(raster::extent(min(alerts_1km$line_x), max(alerts_1km$line_x), min(alerts_1km$line_y), max(alerts_1km$line_y)), "SpatialPolygons");
proj4string(e) = "+proj=longlat +datum=WGS84";

calc_kerneldensity_diggle_512 = function(input,overlay){

	alerts2 = input;

	# Load Waze Alerts data into SpatialPointsDataframe
	alerts2$level=100;

	alerts_level = alerts2[,c("line_x","line_y","level")];
	coordinates(alerts_level) = cbind(alerts_level$line_x , alerts_level$line_y);
	proj4string(alerts_level) = proj4string(e);

	map_wgs84 = spTransform(overlay, CRS("+proj=longlat +datum=WGS84"));

	alerts_level = remove.duplicates(alerts_level);

	window = as.owin(map_wgs84);
	Alerts.ppp = ppp(x=alerts_level@coords[,1],y=alerts_level@coords[,2],window=window,dimyx=c(1024,1024));

	den = density.ppp(Alerts.ppp, sigma = bw.diggle(Alerts.ppp),edge=T);

	den_df = as.data.frame(den);

	return(den);

}

#####################################################
## GENERATE HOTSPOT
#####################################################

for(i in 1:length(subtype)){

	for(j in 1:length(yday)){

		yday_ = yday[j]; 
		subtype_ = subtype[i];

		input = alerts_1km %>% 
		dplyr::select(-Waktu) %>% 
		filter(yday==yday_ & subtype==subtype_) %>%
		as.data.frame();

		output = calc_kerneldensity_diggle_512(input,e);

		nam = paste("hostpot",paste(subtype_, yday_, sep = "_"),sep="_");
 		assign(nam, output);

	}

}

for(i in 1:length(subtype)){

	for(j in 1:length(week)){

		week_ = week[j]; 
		subtype_ = subtype[i];

		input = alerts_1km %>% 
		dplyr::select(-Waktu) %>% 
		filter(week==week_ & subtype==subtype_) %>%
		as.data.frame();

		output = calc_kerneldensity_diggle_512(input,e);

		nam = paste("hostpot",paste(subtype_, week_, sep = "_"),sep="_");
 		assign(nam, output);

	}

}

par(mfrow=c(2,2));

plot(hostpot_ROAD_CLOSED_CONSTRUCTION_40);
plot(hostpot_JAM_STAND_STILL_TRAFFIC_40);
plot(hostpot_HAZARD_WEATHER_FLOOD_40);
plot(hostpot_HAZARD_WEATHER_HAIL_40);

alerts_jams_1km = alerts_1km %>%
dplyr::select(-Waktu) %>%
filter(waze_alerts.type=="JAM") %>%
as.data.frame();

calc_kerneldensity_diggle_512_multi = function(input,overlay){

	alerts2 = input;

	# Load Waze Alerts data into SpatialPointsDataframe
	alerts2$level=100;
	alerts2$level[alerts2$subtype=="JAM_HEAVY_TRAFFIC"]=300;
	alerts2$level[alerts2$subtype=="JAM_STAND_STILL_TRAFFIC"]=500;

	alerts_level = alerts2[,c("line_x","line_y","level")];
	coordinates(alerts_level) = cbind(alerts_level$line_x , alerts_level$line_y);
	proj4string(alerts_level) = proj4string(e);

	map_wgs84 = spTransform(overlay, CRS("+proj=longlat +datum=WGS84"));

	alerts_level = remove.duplicates(alerts_level);

	window = as.owin(map_wgs84);
	Alerts.ppp = ppp(x=alerts_level@coords[,1],y=alerts_level@coords[,2],window=window,dimyx=c(1024,1024));

	den = density.ppp(Alerts.ppp, sigma = bw.diggle(Alerts.ppp),edge=T);

	den_df = as.data.frame(den);

	return(den);

}

	for(j in 1:length(week)){

		week_ = week[j]; 

		input = alerts_1km %>% 
		dplyr::select(-Waktu) %>% 
		filter(week==week_) %>%
		as.data.frame();

		output = calc_kerneldensity_diggle_512_multi(alerts_1km,e);

		nam = paste("hostpot",week_,sep="_");
 		assign(nam, output);

	}

par(mfrow=c(1,4));

plot(hostpot__40);
plot(hostpot__41);
plot(hostpot__42);
plot(hostpot__43);

streets = unique(alerts_1km$waze_alerts.street);

for(i in 1:length(streets)){
	if(i==1){
		output = streets;
	} else{
		output = paste(output,streets,sep="",");
	}
}

#####################################################
## GENERATE HOTSPOT
#####################################################

calc_kerneldensity_diggle_512_multi = function(input,overlay){

	alerts2 = input;

	# Load Waze Alerts data into SpatialPointsDataframe
	alerts2$level=100;
	alerts2$level[alerts2$subtype=="JAM_HEAVY_TRAFFIC"]=300;
	alerts2$level[alerts2$subtype=="JAM_STAND_STILL_TRAFFIC"]=500;

	alerts_level = alerts2[,c("line_x","line_y","level")];
	coordinates(alerts_level) = cbind(alerts_level$line_x , alerts_level$line_y);
	proj4string(alerts_level) = proj4string(e);

	map_wgs84 = spTransform(overlay, CRS("+proj=longlat +datum=WGS84"));

	alerts_level = remove.duplicates(alerts_level);

	window = as.owin(map_wgs84,step=10);
	Alerts.ppp = ppp(x=alerts_level@coords[,1],y=alerts_level@coords[,2],window=window,dimyx=c(512,512));

	den = density.ppp(Alerts.ppp, sigma = bw.diggle(Alerts.ppp),edge=T);

	den_df = as.data.frame(den);

	return(den);

}

output = calc_kerneldensity_diggle_512_multi(alerts_jams_1km,e);

output_df = as.data.frame(output);
den95 = quantile(output_df$value,0.95);
output_95 = output_df %>% filter(value>=den95) %>% as.data.frame();
den95_2 = quantile(output_95$value,0.95);
output_95_2 = output_95 %>% filter(value>=den95_2) %>% as.data.frame();
write.csv(output_df,"jams_pattern_okt.csv",row.names=FALSE);

DBSCAN1 = dbscan(cbind(output_95_2$x, output_95_2$y), eps = 0.00075, minPts = 3);
output_95_2$cluster = DBSCAN1$cluster;

write.csv(output_95_2,"jams_pattern_okt_95.csv",row.names=FALSE);



