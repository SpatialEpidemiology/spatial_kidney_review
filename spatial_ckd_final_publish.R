#installing and loading required packages

listOfPackages <- c("tidyverse","performance","readxl","devtools","gratia"
                    ,"GWmodel","mgcv","RCurl","spgwr","sf","spdep","sp")
for (i in listOfPackages){
  if(! i %in% installed.packages()){
    install.packages(i, dependencies = TRUE)
  }
  require(i)
}

install_github("samclifford/mgcv.helper")

library(tidyverse)
library(readxl)
library(devtools)
library(spgwr)
library(gratia)
library(sf)
library(spdep)
library(sp)
library(mgcv)
library(performance)
library(mgcv.helper)

#read in spatial file
cdh_sf<-st_read("~yourpath/ckd_spatial_final.gpkg")
#converting sf to sp
cdh_sp<-as(cdh_sf,"Spatial")


##creating nearest neighbors matrix via KNN for lagged CKD prev rate
us.nb4<-knearneigh(coordinates(cdh_sp), k=8)
us.nb4<-knn2nb(us.nb4)
us.nb4<-make.sym.nb(us.nb4)
us.wt4<-nb2listw(us.nb4, style="W")


#COUNT MODELING

#poisson model to test for overdispersion
cdh_poi<-mgcv::gam(ckd_ct~scale(dia_rt)+scale(hyp_rt)+scale(EQI)+offset(log(pop65)),family="poisson",data=cdh_sf,method="REML") 
check_overdispersion(cdh_poi)
#overdispersed so use quasi-poisson

#non-spatial quasipoisson GAM
cdh_qp<-gam(ckd_ct~scale(dia_rt)+scale(hyp_rt)+scale(EQI)+offset(log(pop65)),family=quasipoisson(),data=cdh_sf,method="REML")
#vif #no curvature so concurvity() not needed
vif.gam(cdh_qp)
#confidence intervals
con_qp<-mgcv.helper::confint.gam(cdh_qp,levels=0.95)
#confint of prevalence ratios
exp(con_qp$Estimate)
exp(con_qp$`2.5%`)
exp(con_qp$`97.5%`)
#global clustering of residuals
moran.test(residuals(cdh_qp,type="deviance"),listw = us.wt4,alternative="greater")
#adding residuals into dataframe for local moran's I in GeoDa
cdh_sf$cdh_qp_dev_resid<-residuals(cdh_qp,type="deviance")
#informal plotting of residuals
plot(cdh_sf["cdh_qp_dev_resid"],breaks="quantile",border=NA)


#creating lag variable for SLM
cdh_sp$lag_ct<-lag.listw(x=us.wt4, var=cdh_sp$ckd_ct/cdh_sp$pop65)
cdh_sf$lag_ct<-lag.listw(x=us.wt4, var=cdh_sf$ckd_ct/cdh_sf$pop65)
#creating x/y spatial geometry for lag and error count models 
coords<-st_centroid(st_transform(cdh_sf,crs=3857))
coords<-st_transform(coords,crs=4269)
coords<-data.frame(st_coordinates(coords))
cdh_sf$x<-coords$X
cdh_sf$y<-coords$Y


#lagged count model (SLM)
cdh_qp_lag<-gam(ckd_ct~scale(dia_rt)+scale(hyp_rt)+lag_ct+scale(EQI)+offset(log(pop65)),family=quasipoisson(),data=cdh_sf,method="REML")
#vif #no spline curvature so concurvity() not needed
vif.gam(cdh_qp_lag)
#confidence intervals
con_qp_lag<-mgcv.helper::confint.gam(cdh_qp_lag,levels=0.95)
#confint of prevalence ratios #2.5 and 97.5 are just syntax of pkg # they are actually 95%CIs
exp(con_qp_lag$Estimate)
exp(con_qp_lag$`2.5%`)
exp(con_qp_lag$`97.5%`)
#global moran's I
moran.test(residuals(cdh_qp_lag,type="deviance"),listw = us.wt4,alternative="greater")
#adding residuals to dataframe
cdh_sf$cdh_qp_lag_dev_resid<-residuals(cdh_qp_lag,type="deviance")
#informal plotting of residuals
plot(cdh_sf["cdh_qp_lag_dev_resid"],breaks="quantile",border=NA)


#error count model (GAM non-linear spatial terms [thin plate spline tensor]; PSEM)
cdh_qp_te<-gam(ckd_ct~scale(dia_rt)+scale(hyp_rt)+scale(EQI)+offset(log(pop65))+te(x,y,k=14,bs="tp"),family=quasipoisson(),data=cdh_sf,method="REML")
#function to help choose ideal number of spline knots
gam.check(cdh_qp_te)
#vif #concurvity check #no issues detected
vif.gam(cdh_qp_te)
concurvity(cdh_qp_te)
#global moran's I
moran.test(residuals(cdh_qp_te,type="deviance"),listw = us.wt4,alternative="greater")
#adding residuals to dataframe
cdh_sf$cdh_qp_te_dev_resid<-residuals(cdh_qp_te,type="deviance")
#informal plotting of residuals
plot(cdh_sf["cdh_qp_te_dev_resid"],breaks="quantile",border=NA)
#confidence intervals
con_qp_te<-mgcv.helper::confint.gam(cdh_qp_te,levels=0.95)
#confint of prevalence ratios
exp(con_qp_te$Estimate)
exp(con_qp_te$`2.5%`)
exp(con_qp_te$`97.5%`)
#creating tensor product countour map of partial effects
gratia::draw(cdh_qp_te,fun="exp",scales="fixed", n_contour = 8)


#geaographically weighted quasi-Poisson regression
#bandwidth selection #this can take a long time
#using poisson family for bandwidth, as results in same bw as QP
#BANDWIDTH CALCULATION TAKES A LONG WHILE TO RUN
#TO AVOID USE THE OPTION bandwidth=38.7973 in the ggwr() function
#Use same bandwidth in both for better comparison of null w/ adjusted
bw2 <- ggwr.sel(ckd_ct~scale(dia_rt)+scale(EQI)+scale(hyp_rt)+offset(log(pop65)), data=cdh_sp,
                family=poisson(), longlat=FALSE)
#null model for deviance r2 calculation
gwr_null <- ggwr(ckd_ct~1+offset(log(pop65)),
                 data=cdh_sp,  
                 family = quasipoisson(),
                 bandwidth=38.7973,
                 longlat=F,
                 type="deviance")
#main gwqpr model
gwr <- ggwr(ckd_ct~scale(dia_rt)+scale(EQI)+scale(hyp_rt)+offset(log(pop65)),
            data=cdh_sp,  
            family = quasipoisson(),
            bandwidth=38.7973,
            longlat=F,
            type="deviance")
#global moran's I of main gwqpr
moran.test(gwr$SDF$deviance_resids,listw = us.wt4,alternative="greater")
#adding residuals to dataframe
cdh_sf$cdh_qp_gwr_dev_resid<-gwr$SDF$deviance_resids
#informal plotting of residuals
plot(cdh_sf["cdh_qp_gwr_dev_resid"],breaks="quantile",border=NA)

#exp(intercept) and prevalence ratios of gwr output
exp(median(gwr$SDF$X.Intercept.))
exp(min(gwr$SDF$X.Intercept.))
exp(max(gwr$SDF$X.Intercept.))

exp(median(gwr$SDF$scale.dia_rt.))
exp(min(gwr$SDF$scale.dia_rt.))
exp(max(gwr$SDF$scale.dia_rt.))

exp(median(gwr$SDF$scale.EQI.))
exp(min(gwr$SDF$scale.EQI.))
exp(max(gwr$SDF$scale.EQI.))

exp(median(gwr$SDF$scale.hyp_rt.))
exp(min(gwr$SDF$scale.hyp_rt.))
exp(max(gwr$SDF$scale.hyp_rt.))

#adding gwr coefficents to dataframe for mapping
cdh_sf$int_gwr_coef<-exp(gwr$SDF$X.Intercept.)
cdh_sf$EQI_gwr_coef<-exp(gwr$SDF$scale.EQI.)
cdh_sf$dia_gwr_coef<-exp(gwr$SDF$scale.dia_rt.)
cdh_sf$hyp_gwr_coef<-exp(gwr$SDF$scale.hyp_rt.)

#write spatial file with modeling results for external software (GeoDa and ArcGIS)
#shapefile always cuts off variable names
#converting to web mercator projection for better plotting
cdh_sf2<-st_transform(cdh_sf,crs=3857)
#st_write(cdh_sf2,"~yourpath/ckd_spatial_modeling_results.gpkg")
#st_write(cdh_sf2,"~yourpath/ckd_spatial_modeling_results.shp")

#deviance r2 for all models
1-cdh_qp$deviance/cdh_qp$null.deviance
1-cdh_qp_lag$deviance/cdh_qp_lag$null.deviance
1-cdh_qp_te$deviance/cdh_qp_te$null.deviance
1-sum(gwr$SDF$deviance_resids^2)/sum(gwr_null$SDF$deviance_resids^2)

#creating forest plot of prevalence ratios for EQI
dat <- data.frame(
  Index = c(1, 2, 3, 4), ## This provides an order to the data
  label = c("EQI (Non-spatial)", "EQI (SLM)", "EQI (PSEM)", "EQI (GWQPR)"),
  Prevalence_ratio = c(1.013, 1.014, 1.007, 1.012),
  LL = c(1.010, 1.011, 1.004, 1.011),
  UL = c(1.016, 1.016, 1.010, 1.015),
  CI = c("1.010, 1.016", "1.011, 1.016", "1.004, 1.010", "1.011, 1.015")
)

plot1 <- ggplot(dat, aes(y = Index, x = Prevalence_ratio)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = LL, xmax = UL), height = 0.25) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  scale_y_continuous(name = "", breaks=1:4, labels = dat$label, trans = "reverse") +
  xlab("Prevalence ratio 
  (95% CI for QP GAM, SLM, PSEM) 
       (Min-max prevalence ratio for GWQPR)") + 
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x.bottom = element_text(size = 12, colour = "black",vjust = 0.5, hjust=1),
        axis.title.x = element_text(size = 12, colour = "black"))
plot1
