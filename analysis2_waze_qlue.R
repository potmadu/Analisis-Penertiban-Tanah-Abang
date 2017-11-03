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

setwd(working_directory);

#####################################################
## LOAD DATA
#####################################################

qlue = read.csv("qlue_oktober_tanahabang.csv",stringsAsFactors=FALSE);

alerts_1km = read.csv("alerts_1km.csv",stringsAsFactors=FALSE);

alerts_standstill = alerts_1km %>% filter(waze_alerts.subtype=="JAM_STAND_STILL_TRAFFIC") %>% as.data.frame();

alerts_carstopped = read.csv("alerts_carstopped.csv",stringsAsFactors=FALSE);
alerts_pothole = read.csv("alerts_pothole.csv",stringsAsFactors=FALSE);
alerts_construction = read.csv("alerts_construction.csv",stringsAsFactors=FALSE);

alerts1 = alerts_pothole[,c("waze_alerts.uuid","waze_alerts.subtype","line_y","line_x","cluster")];
alerts2 = alerts_construction[,c("waze_alerts.uuid","waze_alerts.subtype","line_y","line_x","cluster")];
alerts3 = alerts_carstopped[,c("waze_alerts.uuid","waze_alerts.subtype","line_y","line_x","cluster")];
alerts4 = alerts_standstill[,c("waze_alerts.uuid","waze_alerts.subtype","line_y","line_x")];

alerts4$cluster=-1;

alerts1$kanal="Waze";
alerts2$kanal="Waze";
alerts3$kanal="Waze";
alerts4$kanal="Waze";

qlue$id_report_source=as.character(qlue$id_report_source);
qlue$cluster=-1;
qlue$kanal="Qlue";

colnames(qlue) = colnames(alerts1);

laporan = rbind(qlue,alerts1,alerts2,alerts3,alerts4);

write.csv(laporan,"waze_qlue.csv",row.names=FALSE);

