library(myr)
library(RNetCDF)

# Load data 
if (FALSE) {
    topo   = my.read.nc("../output/Greenland/GRL-20KM/GRL-20KM_TOPO-RTOPO-2.0.1.nc")

    mar    = my.read.nc("../output/Greenland/GRL-20KM/GRL-20KM_MARv3.5-ERA-30km-ann_1981-2010.nc")
    mar$pr = mar$sf+mar$rf

    rac    = my.read.nc("../output/Greenland/GRL-20KM/GRL-20KM_RACMO23-ERA-INTERIM_ann_1981-2010.nc")
}

if (FALSE) {
    # Check alignment of topopgraphy with our topo dataset 
    image.plot(topo$xc,topo$yc,mar$z_srf-topo$z_srf,zlim=c(-500,500))
    # image.plot(topo$xc,topo$yc,rac$z_srf-topo$z_srf,zlim=c(-100,100))
}

if (FALSE) {
    # Check smb totals to make sure units are ok 
    kk = which(topo$mask >= 2)
    cat("=== Ice sheet annual means [mm/a] (MAR, RACMO) ===","\n")
    cat("Precip  : ", mean(mar$pr[kk]),mean(rac$pr[kk]),"\n")
    cat("Snowfall: ", mean(mar$sf[kk]),mean(rac$sf[kk]),"\n")
    cat("SMB     : ", mean(mar$smb[kk]),mean(rac$smb[kk]),"\n")
}