# This file contains the source code of the functions written during the preparation of
# the Speich et al. article on optimal rooting depth. The purpose of this file is to
# document the steps described in the article and to ensure reproducibility of the results.

# There is NO WARRANTY that this code will be fit for any particular purpose.

# Author: Matthias Speich, WSL

# Contents:
#	- netrad: auxiliary function for the estimation of net radiation
#	- prep.phen: auxiliary function for the calculation of daily phenological status as a function of air temperature
#	- climvars: function to calculate the climate variables required by the optimal rooting depth models
#	- ord.grass: implementation of the optimal rooting depth model of Guswa (2008)
#	- ord: implementation of the optimal rooting depth model of Guswa (2010)
#	- Example application


#############################################

netrad <- function(flx.day,elv){

# Auxiliary function for the estimation of net radiation

# Net radiation, required for the estimation of potential evapotranspiration, is provided in the FLUXNET dataset.
# However, there are many data gaps for this variable. This function estimates net radiation (and its components, 
# net shortwave and net longwave radiation) from air temperature, global radiation and Vapor Pressure Deficit, for which
# temporal coverage is complete in many cases.

# Inputs: 	- flx.day - a data.frame containing daily meteorological measurements, in the format corresponding to the FLUXNET 2015 dataset
#			- elv: elevation in m asl
	
	# Constants
	cp <- 1004
	at <- 0.95
	sigma <- 5.6703E-8
	zeroK <- 273.15
	hundredK <- 373.15
	pZero <- 1013.25
	gammaT <- 0.005
	g <- 9.81
	rAir <- 287.04
	rWat <- 461.5
	rhoWat <- 999.941
	H <- c(13.3185,-1.976,-0.6445,-0.1299)
	pp <- c(-137.0, -75.0,  30.0, 167.0, 236.0, 252.0, 213.0,  69.0, -85.0,-206.0,-256.0,-206.0)
	
	# Select required variables from FLUXNET data.frame
	rg <- flx.day$SW_IN_F			# Global radiation
	rg.pot <- flx.day$SW_IN_POT		# Potential incoming shortwave radiation
	ta <- flx.day$TA_F				# Air temperature
	vpd <- flx.day$VPD_F			# Vapor pressure deficit
	
	
	# Net shortwave radiation
	albedo <- 0.1
	rs <- (1-albedo)*rg
	
	# Auxiliary calculations for net LW: Vapor pressure and relative SSD
	tk <- ta + zeroK
    tr <- 1-(hundredK/tk)
    tp <- tk+gammaT*elv/2
    PT <- pZero * exp(-g*elv/rAir/tp)
    l <- 2500800*(zeroK/tk)
    gamma <- PT*cp/l/(rAir/rWat)
    rhoAir <- 100*PT/rAir/tk
    x<-0
    dx<-0
    for (iter in 1:4){
     	x <- x + H[iter] * tr^iter
    	dx <- dx + iter * H[iter] * tr^(iter-1)
    }
    elx <- pZero * exp(x)
    del <- hundredK*elx/tk/tk*dx
    el <- pmax((elx - vpd), 0.1)
    
    ssd <- pmin(1,(((rg/rg.pot)-0.25)/0.5))
    
    # Net longwave radiation - following Schulla (1997)
    rl <- -sigma*((ta+zeroK)^4.)*(0.52-(0.065*sqrt(el)))* (0.23+(0.77*ssd))  # Schulla 1997
    
    rn <- rs+rl
    
    return(rn)
	
}

#############################################

prep.phen <- function(flx.dd,phen.type,ti=NULL, to=NULL, ta=NULL, cstar=NULL, b=NULL, f=NULL, species, years,allyear=FALSE,months=NULL,outpath,dowrite=FALSE){
	
# Calculates daily phenological status (level of foliage cover for decisuous species, and level of activity for evergreens).
# Based on the model of Kramer (1996) as described by Zierl (2000, 2001) for spring phenology, and on Zierl (2000, 2001) for autumn phenology.

# Inputs:
#	- flx.dd - a data.frame containing daily meteorological measurements, in the format corresponding to the FLUXNET 2015 dataset
#	- phen.type:	1 for deciduous, 2 for evergreen
#	- ti, to, ta, cstar, b, f: parameters of the Kramer model. Not used here, values are attributed based on species (next argument)
# 	- species: allows to specify species name and get the corresponding Kramer parameters
#	- years: vector specifying the years for which phenological status should be calculated
#	- allyear: (logical) - set to TRUE to override calculations and make plants active all year round (i.e. growing season lasts the whole year)
#	- months: vector of month indices (1-12) that can be used to override calculations and specify the months corresponding to the growing season
#	- outpath: path to output file
#	- dowrite: (logical) - set to TRUE to write output to file


# Species-specific parameters:
	if(species=="Picea_abies"){
		ti <- -11.4
		to <- 0.1
		ta <- 16.3
		cstar <- 82.5
		b <- 0.14
		f <- 1.6
		c <- 35.9
	} else if(species=="Fagus_sylvatica"){
		ti <- -21.4
		to <- -1.8
		ta <- 69.7
		cstar <- 115.6
		b <- 0.08
		f <- 2.0
		c <- 47.4
	} else if(species=="Quercus"){
		ti <- -11.4
		to <- -3.8
		ta <- 39.3
		cstar <- 101.7
		b <- 0.11
		f <- 1.9
		c <- 37.8
	} else if(species=="Pinus_sylvestris"){
		ti <- -13.8
		to <- -1.2
		ta <- 16.5
		cstar <- 85.3
		b <- 0.11
		f <- 2.4
		c <- 37.6
	}

	
	flx.dd$year <- as.numeric(substr(flx.day$TIMESTAMP,1,4))
	flx.dd$mon <- as.numeric(substr(flx.day$TIMESTAMP,5,6))
	flx.dd$day <- as.numeric(substr(flx.day$TIMESTAMP,7,8))
	flx.dd <- replace(flx.dd,flx.dd==-9999,NA)
	
	# Subset of the original data.frame to keep only the years of interest (+ one year before, for initialization)
	yrs.df <- c((years[1]-1), years)
	flx.sub <- flx.dd[flx.dd$year %in% yrs.df,]
	
	# Vector with phenological state for each day
	phen.vec <- rep(NA,nrow(flx.sub))
	act.vec <- rep(NA,nrow(flx.sub))
	d1 <- 305  # DOY of November 1st, start of chilling period
	d2 <- 250  # DOY of September 7, start of potential senescence period
	
	# Transform Date to DOY
	datstr <- paste(flx.sub$day,flx.sub$mon,flx.sub$year,sep=".")
	doy <- strftime(as.Date(datstr,format="%d.%m.%Y"),format="%j") # dateToDOY(datstr)
	
	# Start at the first occurrence of Nov 1
	start.days <- which(doy==d1)
	for(i in (1:(length(start.days)-1))){
		# Loop over all years, starting at Nov 1
		
		# State variables for chilling and forcing - initialization at the beginning (Nov 1) of each year
		sc <- 0
		sf <- 0
		
		quies <- FALSE
		active <- FALSE
		senesc <- FALSE
		
		day.bb <- 0
		
		# Indices of start (Nov 1) and end (Oct 31) of each year
		start <- start.days[i]
		end <- (start.days[i+1]-1)
		for(day in start:end){
			
			t <- flx.sub[day,"TA_F"]
			
			if(!quies){
			# Rate of chilling
			    if(t < ti | t > ta){
			    	rc <- 0
			    } else if(t <= to){
				    rc <- (t-ti)/(to-ti)
				} else if(t > to){
					rc <- (t-ta)/(to-ta)
				}
			
				sc <- sc + rc
				if(sc > cstar){quies = TRUE}  # Quiescence phase starts when chilling state reaches a certain threshold c*
				
				phen.vec[day] <- 0
				act.vec[day] <- 0
			}
			else if(!active){
				# During quiescence phase - Rate of forcing is calculated
				if(t > 4.85){
				    rf <- 1/(1+exp(-b*(t-c)))
				 }   else {
				 	rf <- 0
				 }
				sf <- sf + rf
				if(sf > f){
					# Once the state of forcing exceeds the threshold f, 
					active <- TRUE
					day.bb <- day
				}
				
				phen.vec[day] <- 0
				act.vec[day] <- 0
				
			} else if(!senesc) {
				# After bud burst and before onset of leaf senescence
				phen.vec[day] <- min(1,((day-day.bb)/30))
				act.vec[day] <- 1
				# Onset of senescence earliest on DOY 250, after a period of 4 days with mean T < 5°C
				if(doy[day] >= 250){
					mt4days <- mean(flx.sub[(day-3):day,"TA_F"])
					if(mt4days < 5){
						senesc <- TRUE
						day.sen <- day
					}
				}
			} else{
				# Between onset of leaf senescence and start of next chilling period
				phen.vec[day] <- max(0,(1-((day-day.sen)/14)))
				if(phen.vec[day] > 0){
					act.vec[day] <- 1
				} else{
					act.vec[day] <- 0
				}
			}
		
		}

	if(phen.type==2){phen.vec[phen.vec < 1] <- 1}
	
	if(allyear){
		act.vec[1:length(act.vec)] <- 1
	}
	
	# Override all else if months are specified a priori. Used here for Pinus pinaster only, so only dormancy/activity is changed, not leaf phenology
	
	if(!is.null(months)){
		act.vec[1:length(act.vec)] <- 0
		act.vec[flx.sub$mon %in% months] <- 1
	}

	phen.df <- flx.sub[,c("day","mon","year")]
	phen.df$phen <- phen.vec
	phen.df$active <- act.vec
	out.df <- phen.df[phen.df$year %in% years,]
	# Replace NAs with 0 (should occur only for the period after Nov 1 of the last year)
	out.df$phen[is.na(out.df$phen)] <- 0
	out.df$active[is.na(out.df$active)] <- 0
	if(dowrite){write.table(out.df,file=outpath,sep=",",quote=FALSE,row.names=FALSE)}
	}
	
	return(out.df)
}

#############################################

climvars <- function(flx.day,phen,lai,elv,ksimax,kl){
	
# Calculates climate variables required by the optimal rooting depth models of Guswa (2008, 2010)

# Inputs:
#	- flx.day - a data.frame containing daily meteorological measurements, in the format corresponding to the FLUXNET 2015 dataset
#	- phen: a data.frame with information on phenology, each row corresponding to one day. 
#		Must contain a column named "active", with a value greater than 0.0 for days in the growing season, and 
#		equal to 0.0 for days outside the growing season.
#	- LAI: leaf area index during the growing season (numeric)
#	- elv: elevation in m asl (numeric)
#	- ksimax: parameter relating LAI to size of interception reservoir (see Menzel (1997) and Vegas-Galdos et al. (2012))
#	- kl: canopy light extinction coefficient [-]

# Outputs:
#	Returns a list containing the following elements:
#	- pet: mean daily Penman PET for the growing season [mm/day]
#	- prec: mean daily precipitation for the growing season [mm/day]
#	- pfreq: frequency of rainfall events [events/day]
#	- pmdepth: average depth of a rainfall event [mm/event]
#	- ts: mean soil temperature [°C]
#	- mean daily *effective* precipitation (i.e. reaching the ground) for the growing season [mm/day]
# 	- fgs: growing season length, expressed as a fraction of the year
	
	## Set constants
	a <- 86400000
	rAir <- 287.04
	rWat <- 461.5
	rhoWat <- 999.941
	zeroK <- 273.15
	hundredK <- 373.15
	gammaT <- 0.005
	pZero <- 1013.25
	g <- 9.81
	cp <- 1004
	H <- c(13.3185,-1.976,-0.6445,-0.1299)
	pp <- c(-137.0, -75.0,  30.0, 167.0, 236.0, 252.0, 213.0,  69.0, -85.0,-206.0,-256.0,-206.0)
	k <- 0.5
	
	flx.day$year <- as.numeric(substr(flx.day$TIMESTAMP,1,4))
	flx.day$mon <- as.numeric(substr(flx.day$TIMESTAMP,5,6))
	flx.day$day <- as.numeric(substr(flx.day$TIMESTAMP,7,8))
	flx.day <- replace(flx.day,flx.day==-9999,NA)
	
	# Keep only years of selected period (contained in phen file)
	years <- unique(phen$year)
	flx.day <- flx.day[flx.day$year %in% years,]
	
	## Variables from FLUXNET dataset
	t <- flx.day$TA_F				# Air temperature
	tsoil <- flx.day$TS_F_MDS_1		# Soil temperature
	u <- flx.day$WS_F				# Wind speed
	rn <- netrad(flx.day,elv)		# Net radiation (calculated with the auxiliary function netrad, also contained in this file)
	vpd <- flx.day$VPD_F			# Vapor pressure deficit
	
	## Penman potential evaporation
	tk <- t + zeroK
	tr <- 1-(hundredK/tk)
	tp <- tk + gammaT*elv/2
	PT <- pZero * exp(-g*elv/rAir/tp)			# Air pressure
	l <- 2500800*(zeroK/tk)						# Latent heat of vaporization
	gamma <- PT*cp/l/(rAir/rWat)				# Psychrometric "constant"
	rhoAir <- 100*PT/rAir/tk					# Air density
	x<-0
	dx<-0
	for (iter in 1:4){
		x <- x + H[iter] * tr^iter
		dx <- dx + iter * H[iter] * tr^(iter-1)
	}
	elx <- pZero * exp(x)						# Atmospheric water vapor pressure at saturation
	del <- hundredK*elx/tk/tk*dx				# "delta" variable of Penman equation (slope of saturation vapor pressure against temperature)

	# Aerodynamic term of Penman equation, following Penman 1954
	alphaL <- gamma * rhoWat * l/a * 0.263 * (0.5+0.537*u)
	energy <- ((del*rn)+alphaL*vpd)/(gamma+del)
	epen <- energy*a/rhoWat/l					# Penman evaporation
	
	# Calculate sub-canopy fluxes to estimate soil evaporation
	asoil <- pmax(0,rn*(exp(-k*lai)))			# Available energy below the canopy
	pesoil <- ((del*0.8*asoil)/(del+gamma)) *a/rhoWat/l
	
	pet <- vector("numeric",length(years))
	tpot <- vector("numeric",length(years))
	prec <- vector("numeric",length(years))
	pes <- vector("numeric",length(years))
	temp <- vector("numeric",length(years))
	ts <- vector("numeric",length(years))
	for(i in 1:length(years)){
		phen.year <- phen[phen$year==years[i],]
		epen.year <- epen[phen$year==years[i]]
		pet[i] <- mean(epen.year[phen.year$active > 0 & !is.na(phen.year$active)])
		p.year <- flx.day[phen$year==years[i],"P_F"]
		prec[i] <- mean(p.year[phen.year$active > 0 & !is.na(phen.year$active)])
		temp[i] <- mean(t[phen.year$active > 0 & !is.na(phen.year$active)])
		ts[i] <- mean(tsoil[phen.year$active > 0 & !is.na(phen.year$active)])
	}
	
	## Temporal distribution of precipitation
	pmdepth <- vector("numeric",length(years))
	pfreq <- vector("numeric",length(years))
	for(i in 1:length(years)){
	# Get precipitation of the current growing season
	gsprec <- flx.day[phen$year==years[i] & phen$active > 0, "P_F"]
	prevp <- FALSE
	nevts <- 0
	 for(j in 1:length(gsprec)){
	 	if(gsprec[j] >= 0.5 & !is.na(gsprec[j])){
			if(!prevp){nevts <- nevts+1}
			prevp <- TRUE
		} else{
			prevp <- FALSE
	 	}
	 }
	 pfreq[i] <- nevts/length(gsprec)
	 pmdepth[i] <- sum(gsprec/nevts)
	}
	
	climate.df <- data.frame(years,pet,prec,pfreq,pmdepth,ts)
	clim.mean <- apply(climate.df,2,mean)	# Get mean values over all available years
	
	# Effective precipitation - see Guswa (2008), Eq. 6
	# Assumption: interception evaporative depth corresponds to canopy interception capacity
	simax <- ksimax*log10(1+lai)
	clim.mean$peff <- max(clim.mean[3] * exp(-simax/clim.mean[5]))
	
	# Growing season length (fraction of a year)
	clim.mean$fgs <- sum(phen$active)/nrow(phen)
	
	# Potential transpiration 
	clim.mean$tpot <- clim.mean$pet * (1-exp(-kl*lai)) * 0.75
	
	# Potential soil/understory transpiration
	clim.mean$pes <- clim.mean$pet * exp(-kl*lai)
	
	clim.mean <- clim.mean[2:length(clim.mean)]
	
	return(clim.mean)
}

#############################################

ord.grass <- function(w,tpot,prec,avdepth,whc,rresp,rld,srl,wue,fgs){
	
# Implementation of the optimal rooting depth function of Guswa (2008)

# Inputs:
# 	- tpot: mean daily potential transpiration in the growing season [mm/day]
# 	- prec: mean daily precipitation in the growing season [mm/day]
# 	- avdepth: Average depth of a precipitation event [mm/day]
# 	- whc: Soil water holding capacity [mm/mm]
# 	- rresp: Root respiration rate
# 	- rld: Root length density
# 	- srl: specific rooting length
# 	- wue: Water-use efficiency
# 	- fgs: fraction of the year occupied by GS

# Outputs: 
#	Returns a vector containing 3 values:
#	1: optimal rooting zone storage capacity, normalized by the average depth of a rainfall event (see e.g. Porporato et al. 2004)
#	2: optimal rooting depth [mm], as per Guswa (2008)
#	3: optimal rooting zone storage capacity [mm, i.e. l/m2]
	

	
	out <- vector("numeric",3)
	
	a <- (rresp*rld)/(srl*wue*tpot*fgs)
	
	if(w>=1){
		x <- w*((1+(whc/avdepth)*(((1-w)^2)/(2*a))) - sqrt((whc/avdepth)*(((1-w)^2)/(a)) + (((whc/avdepth)*(((1-w)^2)/(2*a)))^2)))
	} else {
		x <- w*((1+(whc/avdepth)*(((1-w)^2)/(2*a))) + sqrt((whc/avdepth)*(((1-w)^2)/(a)) + (((whc/avdepth)*(((1-w)^2)/(2*a)))^2)))
	}
	
	zr <- (avdepth/(whc*(1-w)))*log(x)
	
	out[1] <- zr/(avdepth/whc)
	out[2] <- zr
	out[3] <- zr*whc
	
	return(out)
}

#############################################


ord <- function(w,tpot,prec,avdepth,whc,fgs,wue,gr,srl,rld){
	
# New implementation of Guswa's 2010 model
# Older version (provided as supplement of the HESS discussion paper) gave incorrect results

# Inputs:
# 	- tpot: mean daily potential transpiration in the growing season [mm/day]
# 	- prec: mean daily precipitation in the growing season [mm/day]
# 	- avdepth: Average depth of a precipitation event [mm/day]
# 	- whc: Soil water holding capacity [mm/mm]
# 	- rresp: Root respiration rate
# 	- rld: Root length density
# 	- srl: specific rooting length
# 	- wue: Water-use efficiency
# 	- fgs: fraction of the year occupied by GS

# Outputs: 
#	Returns a vector containing 3 values:
#	1: optimal rooting zone storage capacity, normalized by the average depth of a rainfall event (see e.g. Porporato et al. 2004)
#	2: optimal rooting depth [mm], as per Guswa (2008)
#	3: optimal rooting zone storage capacity [mm, i.e. l/m2]
	
	out <- vector("numeric",3)
		
	para <- (gr*rld)/(srl * wue * fgs)
	
	xvals <- seq(0,3000,by=0.1)
	porp.pred <- porp.ze(xvals, avdepth, whc, w, tpot)
	# Ze: effective rooting depth
	
	porp.diff <- diff(porp.pred) * 10
	
	z <- xvals[which(porp.diff <= para)[1]]
	
	out[1] <- z/(avdepth/whc)
	out[2] <- z
	out[3] <- z * whc
	
	return(out)
	
}

#############################################

# Example application


### Constant vegetation parameters (following Donohue et al. 2012)

# Grass parameters
gr20.grass <- 0.5
rld.grass <- 0.1
srl.grass <- 1500
wue.grass <- 0.22
fgs.grass <- 0.7
kl <- 0.5

# Tree parameters
gr20 <- 0.5
rld <- 0.1
srl <- 1500
wue <- 0.33



# Site parameters (corresponding to the site Tharandt)
lai <- 7.2
whc <- 0.154
elv <- 320
	


# Read in a FLUXNET daily file using read.table() as a data.frame named flx.day
# flx.day <- read.table(...)

# Calculate phenological status
#phen <- prep.phen(flx.dd=flx.day,phen.type=2,species="Pinus_sylvestris",years=2003:2007,outpath="dummy")
phen <- prep.phen(flx.dd=flx.day,phen.type=2,species="Picea_abies",years=1997:2003,outpath="dummy")

# Get climate variables (warning is normal)
clim <- climvars(flx.day,phen,lai=lai,elv=elv,ksimax=2,kl=0.5) 

### Calculate optimal rooting depth for understory and canopy
pet <- clim$pet
prec <- clim$peff

pmdepth <- clim$pmdepth
fgs <- clim$fgs
ts <- clim$ts

# Root respiration rate as a function of temperature
gr <- gr20*2^((ts-20)/10)
gr.grass <- gr20.grass*2^((ts-20)/10)

# Data.frame containing results
sfc.df <- data.frame(zno=numeric(),zo=numeric(),sfco=numeric(),znu=numeric(),zu=numeric(),sfcu=numeric())

sfc.df[1,1:3] <- ord(w,(pet*0.75),prec,pmdepth,whc,fgs,wue,gr,srl,rld)
sfc.df[1,4:6] <- ord.grass(w,(pet*0.75),prec,pmdepth,whc,gr.grass,rld.grass,srl.grass,wue.grass,fgs.grass)

# Add understory and overstory Sr to get the total Sr:

sr <- sfc.df$sfco * (1-exp(-kl*lai))  + sfc.df$sfcu * exp(-kl*lai)

