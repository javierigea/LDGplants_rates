##JI: this is a modified version of make_models_estimate_noroot.R from Caetano et al Evolution 2018

## This is the set of models we are using to estimate.

## Note that we have a block of GeoSSE models, then fGeoHiSSE models and then MuSSE models.
## The models are set so they match in parameter complexity.

#run model in hydrogen cluster
.libPaths('/home/ji247/R/x86_64-pc-linux-gnu-library/3.5')
library(hisse)

## The function to run the models.

geohisse.ldgplants.evaluate.models <- function(model.number){
    ## Function to evaluate all the models in the set and save the output.
    ## We will check all the models for each replicate in sequence and each replicate in parallel.
    ## phy = phylo
    ## dat = matrix with the data
    ## outname = a string with the name for the replicate.
    GBIFdata<-read.csv('/home/ji247/ldg_plants/raw_data/GBIFdatasummary.csv')
	# For GeoHiSSE.First column has the species names and second column has area codes. Values for the areas need to be 0, 1, or 2, where 0 is the widespread area '01', 1 is endemic area '00' and 2 is endemic area '11'. See 'Details'.
	#strict tropical 2 = abs(max latitude)<23.5 & abs(min latitude)<23.5  & median.latitude <23.5
	GBIFdata[abs(GBIFdata$Max.Latitude)<=23.5&abs(GBIFdata$Min.Latitude)<=23.5&abs(GBIFdata$Median.Latitude)<=23.5,'strict.tropical']<-2
	#strict tropical 1 (strict temperate)
	GBIFdata[abs(GBIFdata$Max.Latitude)>23.5&abs(GBIFdata$Min.Latitude)>23.5&abs(GBIFdata$Median.Latitude)>23.5,'strict.tropical']<-1
	#0: not strict tropical or strict temperate species
	GBIFdata$strict.tropical[is.na(GBIFdata$strict.tropical)]<-0
	GBIFdata$binomial<-paste(GBIFdata$Genus.Name,GBIFdata$Species.Name,sep='_')
	table(GBIFdata$strict.tropical)
	
	GBOTB.bif.angios_BAMM.tree<-read.tree('/home/ji247/ldg_plants/GBOTB.bif.angios_BAMM.tree')
	GBIFdata.BAMM<-GBIFdata[GBIFdata$binomial%in%GBOTB.bif.angios_BAMM.tree$tip.label,]
	GBIFdata.BAMM<-GBIFdata.BAMM[-which(duplicated(GBIFdata.BAMM$binomial)),]
	GBOBT_GBIFdata<-drop.tip(GBOTB.bif.angios_BAMM.tree,setdiff(GBOTB.bif.angios_BAMM.tree$tip.label,GBIFdata.BAMM$binomial))
	traits<-GBIFdata.BAMM$strict.tropical
	names(traits)<-GBIFdata.BAMM$binomial
	traits.df<-as.data.frame(cbind(names(traits),unname(traits)),stringsAsFactors=F)
	colnames(traits.df)<-c('species','trait')
	sampling0<-0.2417261
	sampling1<-0.3013663
	sampling2<-0.1258426
	phy<-GBOBT_GBIFdata
	dat<-traits.df
    
    ###############################################################################
    ## Block of GeoSSE-like models.
    ## Here extirpation is linked to range reduction (change from AB to A or to B).
    ## Also jumps are not allowed.
    ###############################################################################
    
    ## Model 1 - Dispersal parameters vary only, no range-dependent diversification.
    if(model.number == 1){
        speciation <- c(1,1,1)
        extirpation <- c(1,1)
        trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=0, include.jumps=FALSE, separate.extirpation=FALSE)
        mod1 <- try( fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
        , hidden.areas=FALSE, trans.rate=trans.rate, assume.cladogenetic=TRUE) )
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        saveRDS(mod1, file=paste0('ldg_geohisse', "_mod1.rds"))
        
    }
    ## Model 2. Canonical GeoSSE model, range effect on diversification
    if(model.number == 2){
        speciation <- c(1,2,3)
        extirpation <- c(1,2)
        trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=0, include.jumps=FALSE, separate.extirpation=FALSE)
        mod2 <- try( fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
        , hidden.areas=FALSE, trans.rate=trans.rate, assume.cladogenetic=TRUE) )
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        saveRDS(mod2, file=paste0('ldg_geohisse', "_mod2.rds"))
        
    }
    ## Model 3. Heterogeneous diversification, not tied to range evolution.
    ## Assumes three distinct diversification rates.
    if(model.number == 3){
        speciation <- c(1,1,1,2,2,2,3,3,3)
        extirpation <- c(1,1,2,2,3,3)
        trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=2, make.null=TRUE, include.jumps=FALSE, separate.extirpation=FALSE)
        mod3 <- try( fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
        , hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE) )
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        saveRDS(mod3, file=paste0('ldg_geohisse', "_mod3.rds"))
        

    }
    ## Model 4. Heterogeneous diversification, tied to range evolution.
    ## Assumes 6 distinct diversification rates.
    if(model.number == 4){
    	speciation <- c(1,2,3,4,5,6)
    	extirpation <- c(1,2,3,4)
    	trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=1, include.jumps=FALSE, separate.extirpation=FALSE)
    	mod4 <- try( fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
    	                 , hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE) )
    	capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    	saveRDS(mod4, file=paste0('ldg_geohisse', "_mod4.rds"))

    }
    ## Model 5. Heterogeneous diversification, not tied to range evolution. Assumes 5 distinct diversification rates.
    if(model.number == 5){
        speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3))
        extirpation <- c(rep(1,2), rep(2,2), rep(3,2), rep(4,2), rep(5,2))
        trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=4, make.null=TRUE, include.jumps=FALSE, separate.extirpation=FALSE)
        mod5 <- try( fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
        , hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE) )
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        saveRDS(mod5, file=paste0('ldg_geohisse', "_mod5.rds"))
    }
    ## Model 6. Heterogeneous diversification, not tied to range evolution. Assumes two distinct diversification rates.
    if(model.number == 6){
        speciation <- c(1,1,1,2,2,2)
        extirpation <- c(1,1,2,2)
        trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=1, make.null=TRUE, include.jumps=FALSE, separate.extirpation=FALSE)
        mod6 <- try( fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
        , hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE) )
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        saveRDS(mod6, file=paste0('ldg_geohisse', "_mod6.rds"))
    }
    ###############################################################################
    ## Block of cladogenetic models not GeoSSE-like.
    ## Here extirpation is NOT linked to range reduction.
    ## So range reduction is different from the extinction of an endemic lineage.
    ## Jumps between endemic areas are not allowed.
    ###############################################################################
    
    ## Model 1 - Dispersal parameters vary only, no range-dependent diversification.
    if(model.number == 7){
        speciation <- c(1,1,1)
        extirpation <- c(1,1)
        trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=0, include.jumps=FALSE, separate.extirpation=TRUE)
        mod7 <- try( fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
        , hidden.areas=FALSE, trans.rate=trans.rate, assume.cladogenetic=TRUE) )
        saveRDS(mod7, file=paste0('ldg_geohisse', "_mod7.rds"))
    }
    ## Model 2. Canonical GeoSSE model, range effect on diversification
    if(model.number == 8){
        print(8)
        speciation <- c(1,2,3)
        extirpation <- c(1,2)
        trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=0, include.jumps=FALSE, separate.extirpation=TRUE)
        mod8 <- try( fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
        , hidden.areas=FALSE, trans.rate=trans.rate, assume.cladogenetic=TRUE) )
        saveRDS(mod8, file=paste0('ldg_geohisse', "_mod8.rds"))
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        capture.output( print( paste0("Analysis - ", 'ldg_geohisse', " - model 8 done.") ), file=paste0('ldg_geohisse',".log"), append=TRUE )
    }
    ## Model 3. Heterogeneous diversification, not tied to range evolution.
    ## Assumes three distinct diversification rates.
    if(model.number == 9){
        speciation <- c(1,1,1,2,2,2,3,3,3)
        extirpation <- c(1,1,2,2,3,3)
        trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=2, make.null=TRUE, include.jumps=FALSE, separate.extirpation=TRUE)
        mod9 <- try( fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
        , hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE) )
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        saveRDS(mod9, file=paste0('ldg_geohisse', "_mod9.rds"))
    }
    ## Model 4. Heterogeneous diversification, tied to range evolution.
    ## Assumes 6 distinct diversification rates.
    if(model.number == 10){
        speciation <- c(1,2,3,4,5,6)
        extirpation <- c(1,2,3,4)
        trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=1, include.jumps=FALSE, separate.extirpation=TRUE)
        mod10 <- try( fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
        , hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE) )
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        saveRDS(mod10, file=paste0('ldg_geohisse', "_mod10.rds"))
    }
    ## Model 5. Heterogeneous diversification, not tied to range evolution. Assumes 5 distinct diversification rates.
    if(model.number == 11){
        speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3))
        extirpation <- c(rep(1,2), rep(2,2), rep(3,2), rep(4,2), rep(5,2))
        trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=4, make.null=TRUE, include.jumps=FALSE, separate.extirpation=TRUE)
        mod11 <- try( fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
        , hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE) )
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        saveRDS(mod11, file=paste0('ldg_geohisse', "_mod11.rds"))
    }
    ## Model 6. Heterogeneous diversification, not tied to range evolution. Assumes two distinct diversification rates.
    if(model.number == 12){
        speciation <- c(1,1,1,2,2,2)
        extirpation <- c(1,1,2,2)
        trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=1, make.null=TRUE, include.jumps=FALSE, separate.extirpation=TRUE)
        mod12 <- try( fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
        , hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE) )
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        saveRDS(mod12, file=paste0('ldg_geohisse', "_mod12.rds"))
    }
    
    ###############################################################################
    ## Second block of anagenetic models MuSSE-like.
    ## These are very liberal models. They are really GeoSSE-like without cladogenetic events
    ###############################################################################

    ## Model 1. Transitions only. No character effect on diversification
    if(model.number == 13){
        speciation <- c(1,1,1)
        extirpation <- c(1,1,1)
        trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=0, include.jumps=FALSE, separate.extirpation=TRUE)
        mod13 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
        , hidden.areas=FALSE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        saveRDS(mod13, file=paste0('ldg_geohisse', "_mod13.rds"))
    
    }
    ## Model 2. Character effect on diversification.
    if(model.number == 14){
        speciation <- c(1,2,3)
        extirpation <- c(1,2,3)
        trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=0, include.jumps=FALSE, separate.extirpation=TRUE)
        mod14 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
        , hidden.areas=FALSE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        saveRDS(mod14, file=paste0('ldg_geohisse', "_mod14.rds"))
   
    }
    ## Model 3. No character effect on diversification.
    if(model.number == 15){
        speciation <- c(1,1,1,2,2,2,3,3,3)
        extirpation <- c(1,1,1,2,2,2,3,3,3)
        trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=2, include.jumps=FALSE, separate.extirpation=TRUE, make.null=TRUE)
        mod15 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
        , hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        saveRDS(mod15, file=paste0('ldg_geohisse', "_mod15.rds"))
    
    }
    ## Model 4. Character effect on diversification, with a hidden state
    if(model.number == 16){
        speciation <- c(1,2,3,4,5,6)
        extirpation <- c(1,2,3,4,5,6)
        trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=1, include.jumps=FALSE, separate.extirpation=TRUE)
        mod16 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation, hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        saveRDS(mod16, file=paste0('ldg_geohisse', "_mod16.rds"))
   
    }
    ## Model 5. No character effect on diversification, multiple shifts
    if(model.number == 17){
        speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3))
        extirpation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3))
        trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=4, include.jumps=FALSE, separate.extirpation=TRUE, make.null=TRUE)
        mod17 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation, hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        saveRDS(mod17, file=paste0('ldg_geohisse', "_mod17.rds"))
   
    }
    ## Model 6*. No character effect on diversification, multiple shifts.
    if(model.number == 18){
        speciation <- c(rep(1,3), rep(2,3))
        extirpation <- c(rep(1,3), rep(2,3))
        trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=1, include.jumps=FALSE, separate.extirpation=TRUE, make.null=TRUE)
        mod18 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation, hidden.areas=TRUE
        , trans.rate=trans.rate, assume.cladogenetic=FALSE)
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        saveRDS(mod18, file=paste0('ldg_geohisse', "_mod18.rds"))
  
    }
    ## Just trying to complete the sheeit -- four character independent shifts
    if(model.number == 19){
        speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3))
        extirpation <- c(rep(1,2), rep(2,2), rep(3,2), rep(4,2))
        trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=3, include.jumps=FALSE, separate.extirpation=FALSE, make.null=TRUE)
        mod19 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation, hidden.areas=TRUE
        , trans.rate=trans.rate, assume.cladogenetic=TRUE)
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        saveRDS(mod19, file=paste0('ldg_geohisse', "_mod19.rds"))

    }
    if(model.number == 20){
        speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3))
        extirpation <- c(rep(1,2), rep(2,2), rep(3,2), rep(4,2))
        trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=3, include.jumps=FALSE, separate.extirpation=TRUE, make.null=TRUE)
        mod20 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation, hidden.areas=TRUE
        , trans.rate=trans.rate, assume.cladogenetic=TRUE)
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        saveRDS(mod20, file=paste0('ldg_geohisse', "_mod20.rds"))

    }
    if(model.number == 21){
        speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3))
        extirpation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3))
        trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=3, include.jumps=FALSE, separate.extirpation=TRUE, make.null=TRUE)
        mod21 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation, hidden.areas=TRUE
        , trans.rate=trans.rate, assume.cladogenetic=FALSE)
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        saveRDS(mod21, file=paste0('ldg_geohisse', "_mod21.rds"))

    }
    ## A special type of MuSSE model: No speciation or extinction in the widespread range. These are endemic only linked processes. We assume xA and xB are linked to transitions from AB.
    ## Model 1. Transitions only. No character effect on diversification
    if(model.number == 22){
        speciation <- c(1,1,0)
        extirpation <- c(1,1,0)
        trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=0, include.jumps=FALSE, separate.extirpation=FALSE)
        mod22 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
        , hidden.areas=FALSE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        saveRDS(mod22, file=paste0('ldg_geohisse', "_mod22.rds"))
 
    }
    ## Model 2. Character effect on diversification.
    if(model.number == 23){
        speciation <- c(1,2,0)
        extirpation <- c(1,2,0)
        trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=0, include.jumps=FALSE, separate.extirpation=FALSE)
        mod23 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
        , hidden.areas=FALSE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
        saveRDS(mod23, file=paste0('ldg_geohisse', "_mod23.rds"))
 
    }
    ## Model 3. No character effect on diversification.
    if(model.number == 24){
        speciation <- c(1,1,0,2,2,0,3,3,0)
        extirpation <- c(1,1,0,2,2,0,3,3,0)
        trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=2, include.jumps=FALSE, separate.extirpation=FALSE, make.null=TRUE)
        mod24 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
        , hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        saveRDS(mod24, file=paste0('ldg_geohisse', "_mod24.rds"))

    }
    ## Model 4. Character effect on diversification, with a hidden state
    if(model.number == 25){
        speciation <- c(1,2,0,3,4,0)
        extirpation <- c(1,2,0,3,4,0)
        trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=1, include.jumps=FALSE, separate.extirpation=FALSE)
        mod25 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation, hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
        saveRDS(mod25, file=paste0('ldg_geohisse', "_mod25.rds"))

    }
    ## Model 5. No character effect on diversification, multiple shifts
    if(model.number == 26){
        speciation <- c(rep(1,2),0, rep(2,2),0, rep(3,2),0, rep(4,2),0, rep(5,2),0)
        extirpation <- c(rep(1,2),0, rep(2,2),0, rep(3,2),0, rep(4,2),0, rep(5,2),0)
        trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=4, include.jumps=FALSE, separate.extirpation=FALSE, make.null=TRUE)
        mod26 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation, hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        saveRDS(mod26, file=paste0('ldg_geohisse', "_mod26.rds"))
 
    }
    ## Model 6*. No character effect on diversification, multiple shifts.
    if(model.number == 27){
        speciation <- c(rep(1,2),0, rep(2,2),0)
        extirpation <- c(rep(1,2),0, rep(2,2),0)
        trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=1, include.jumps=FALSE, separate.extirpation=FALSE, make.null=TRUE)
        mod27 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation, hidden.areas=TRUE
        , trans.rate=trans.rate, assume.cladogenetic=FALSE)
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        saveRDS(mod27, file=paste0('ldg_geohisse', "_mod27.rds"))
 
    }
    if(model.number == 28){
        speciation <- c(rep(1,2),0, rep(2,2),0, rep(3,2),0, rep(4,2),0)
        extirpation <- c(rep(1,2),0, rep(2,2),0, rep(3,2),0, rep(4,2),0)
        trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=3, include.jumps=FALSE, separate.extirpation=FALSE, make.null=TRUE)
        mod28 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation, hidden.areas=TRUE
        , trans.rate=trans.rate, assume.cladogenetic=FALSE)
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        saveRDS(mod28, file=paste0('ldg_geohisse', "_mod28.rds"))
 
    }
   ## A special type of MuSSE model: No speciation or extinction in the widespread range. These are endemic only linked processes. We assume xA and xB are unlinked to transitions from AB.
    if(model.number == 29){
        speciation <- c(1,1,0)
        extirpation <- c(1,1,0)
        trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=0, include.jumps=FALSE, separate.extirpation=TRUE)
        mod29 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
        , hidden.areas=FALSE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        saveRDS(mod29, file=paste0('ldg_geohisse', "_mod29.rds"))
 
    }
    ## Model 2. Character effect on diversification.
    if(model.number == 30){
        speciation <- c(1,2,0)
        extirpation <- c(1,2,0)
        trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=0, include.jumps=FALSE, separate.extirpation=TRUE)
        mod30 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
        , hidden.areas=FALSE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        saveRDS(mod30, file=paste0('ldg_geohisse', "_mod30.rds"))
   
    }
    ## Model 3. No character effect on diversification.
    if(model.number == 31){
        speciation <- c(1,1,0,2,2,0,3,3,0)
        extirpation <- c(1,1,0,2,2,0,3,3,0)
        trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=2, include.jumps=FALSE, separate.extirpation=TRUE, make.null=TRUE)
        mod31 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
        , hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
        saveRDS(mod31, file=paste0('ldg_geohisse', "_mod31.rds"))
    
    }
    ## Model 4. Character effect on diversification, with a hidden state
    if(model.number == 32){
        speciation <- c(1,2,0,3,4,0)
        extirpation <- c(1,2,0,3,4,0)
        trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=1, include.jumps=FALSE, separate.extirpation=TRUE)
        mod32 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation, hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        saveRDS(mod32, file=paste0('ldg_geohisse', "_mod32.rds"))
   
    }
    ## Model 5. No character effect on diversification, multiple shifts
    if(model.number == 33){
        speciation <- c(rep(1,2),0, rep(2,2),0, rep(3,2),0, rep(4,2),0, rep(5,2),0)
        extirpation <- c(rep(1,2),0, rep(2,2),0, rep(3,2),0, rep(4,2),0, rep(5,2),0)
        trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=4, include.jumps=FALSE, separate.extirpation=TRUE, make.null=TRUE)
        mod33 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation, hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        saveRDS(mod33, file=paste0('ldg_geohisse', "_mod33.rds"))
  
    }
    ## Model 6*. No character effect on diversification, multiple shifts.
    if(model.number == 34){
        speciation <- c(rep(1,2),0, rep(2,2),0)
        extirpation <- c(rep(1,2),0, rep(2,2),0)
        trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=1, include.jumps=FALSE, separate.extirpation=TRUE, make.null=TRUE)
        mod34 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation, hidden.areas=TRUE
        , trans.rate=trans.rate, assume.cladogenetic=FALSE)
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        saveRDS(mod34, file=paste0('ldg_geohisse', "_mod34.rds"))
   
    }
    if(model.number == 35){
        speciation <- c(rep(1,2),0, rep(2,2),0, rep(3,2),0, rep(4,2),0)
        extirpation <- c(rep(1,2),0, rep(2,2),0, rep(3,2),0, rep(4,2),0)
        trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=3, include.jumps=FALSE, separate.extirpation=TRUE, make.null=TRUE)
        mod35 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation, hidden.areas=TRUE
        , trans.rate=trans.rate, assume.cladogenetic=FALSE)
        capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
        saveRDS(mod35, file=paste0('ldg_geohisse', "_mod35.rds"))
  
    }
    

}

geohisse.ldgplants.evaluate.models_local <- function(model.number){
  ## Function to evaluate all the models in the set and save the output.
  ## We will check all the models for each replicate in sequence and each replicate in parallel.
  ## phy = phylo
  ## dat = matrix with the data
  ## outname = a string with the name for the replicate.
  GBIFdata<-read.csv('~/Dropbox/Work_in_progress/LDG_plants/raw_data/GBIFdatasummary.csv')
  # For GeoHiSSE.First column has the species names and second column has area codes. Values for the areas need to be 0, 1, or 2, where 0 is the widespread area '01', 1 is endemic area '00' and 2 is endemic area '11'. See 'Details'.
  #strict tropical 2 = abs(max latitude)<23.5 & abs(min latitude)<23.5  & median.latitude <23.5
  GBIFdata[abs(GBIFdata$Max.Latitude)<=23.5&abs(GBIFdata$Min.Latitude)<=23.5&abs(GBIFdata$Median.Latitude)<=23.5,'strict.tropical']<-2
  #strict tropical 1 (strict temperate)
  GBIFdata[abs(GBIFdata$Max.Latitude)>23.5&abs(GBIFdata$Min.Latitude)>23.5&abs(GBIFdata$Median.Latitude)>23.5,'strict.tropical']<-1
  #0: not strict tropical or strict temperate species
  GBIFdata$strict.tropical[is.na(GBIFdata$strict.tropical)]<-0
  GBIFdata$binomial<-paste(GBIFdata$Genus.Name,GBIFdata$Species.Name,sep='_')
  table(GBIFdata$strict.tropical)
  
  GBOTB.bif.angios_BAMM.tree<-read.tree('~/Dropbox/Work_in_progress/LDG_plants/raw_data/GBOTB.bif.angios_BAMM.tree')
  GBIFdata.BAMM<-GBIFdata[GBIFdata$binomial%in%GBOTB.bif.angios_BAMM.tree$tip.label,]
  GBIFdata.BAMM<-GBIFdata.BAMM[-which(duplicated(GBIFdata.BAMM$binomial)),]
  GBOBT_GBIFdata<-drop.tip(GBOTB.bif.angios_BAMM.tree,setdiff(GBOTB.bif.angios_BAMM.tree$tip.label,GBIFdata.BAMM$binomial))
  traits<-GBIFdata.BAMM$strict.tropical
  names(traits)<-GBIFdata.BAMM$binomial
  traits.df<-as.data.frame(cbind(names(traits),unname(traits)),stringsAsFactors=F)
  colnames(traits.df)<-c('species','trait')
  sampling0<-0.2417261
  sampling1<-0.3013663
  sampling2<-0.1258426
  phy<-GBOBT_GBIFdata
  dat<-traits.df
  
  ###############################################################################
  ## Block of GeoSSE-like models.
  ## Here extirpation is linked to range reduction (change from AB to A or to B).
  ## Also jumps are not allowed.
  ###############################################################################
  
  ## Model 1 - Dispersal parameters vary only, no range-dependent diversification.
  if(model.number == 1){
    speciation <- c(1,1,1)
    extirpation <- c(1,1)
    trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=0, include.jumps=FALSE, separate.extirpation=FALSE)
    mod1 <- try( fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
                           , hidden.areas=FALSE, trans.rate=trans.rate, assume.cladogenetic=TRUE) )
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod1, file=paste0('ldg_geohisse', "_mod1.rds"))
    
  }
  ## Model 2. Canonical GeoSSE model, range effect on diversification
  if(model.number == 2){
    speciation <- c(1,2,3)
    extirpation <- c(1,2)
    trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=0, include.jumps=FALSE, separate.extirpation=FALSE)
    mod2 <- try( fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
                           , hidden.areas=FALSE, trans.rate=trans.rate, assume.cladogenetic=TRUE) )
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod2, file=paste0('ldg_geohisse', "_mod2.rds"))
    
  }
  ## Model 3. Heterogeneous diversification, not tied to range evolution.
  ## Assumes three distinct diversification rates.
  if(model.number == 3){
    speciation <- c(1,1,1,2,2,2,3,3,3)
    extirpation <- c(1,1,2,2,3,3)
    trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=2, make.null=TRUE, include.jumps=FALSE, separate.extirpation=FALSE)
    mod3 <- try( fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
                           , hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE) )
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod3, file=paste0('ldg_geohisse', "_mod3.rds"))
    
    
  }
  ## Model 4. Heterogeneous diversification, tied to range evolution.
  ## Assumes 6 distinct diversification rates.
  if(model.number == 4){
    speciation <- c(1,2,3,4,5,6)
    extirpation <- c(1,2,3,4)
    trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=1, include.jumps=FALSE, separate.extirpation=FALSE)
    mod4 <- try( fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
                           , hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE) )
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod4, file=paste0('ldg_geohisse', "_mod4.rds"))
    
  }
  ## Model 5. Heterogeneous diversification, not tied to range evolution. Assumes 5 distinct diversification rates.
  if(model.number == 5){
    speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3))
    extirpation <- c(rep(1,2), rep(2,2), rep(3,2), rep(4,2), rep(5,2))
    trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=4, make.null=TRUE, include.jumps=FALSE, separate.extirpation=FALSE)
    mod5 <- try( fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
                           , hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE) )
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod5, file=paste0('ldg_geohisse', "_mod5.rds"))
  }
  ## Model 6. Heterogeneous diversification, not tied to range evolution. Assumes two distinct diversification rates.
  if(model.number == 6){
    speciation <- c(1,1,1,2,2,2)
    extirpation <- c(1,1,2,2)
    trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=1, make.null=TRUE, include.jumps=FALSE, separate.extirpation=FALSE)
    mod6 <- try( fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
                           , hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE) )
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod6, file=paste0('ldg_geohisse', "_mod6.rds"))
  }
  ###############################################################################
  ## Block of cladogenetic models not GeoSSE-like.
  ## Here extirpation is NOT linked to range reduction.
  ## So range reduction is different from the extinction of an endemic lineage.
  ## Jumps between endemic areas are not allowed.
  ###############################################################################
  
  ## Model 1 - Dispersal parameters vary only, no range-dependent diversification.
  if(model.number == 7){
    speciation <- c(1,1,1)
    extirpation <- c(1,1)
    trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=0, include.jumps=FALSE, separate.extirpation=TRUE)
    mod7 <- try( fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
                           , hidden.areas=FALSE, trans.rate=trans.rate, assume.cladogenetic=TRUE) )
    saveRDS(mod7, file=paste0('ldg_geohisse', "_mod7.rds"))
  }
  ## Model 2. Canonical GeoSSE model, range effect on diversification
  if(model.number == 8){
    print(8)
    speciation <- c(1,2,3)
    extirpation <- c(1,2)
    trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=0, include.jumps=FALSE, separate.extirpation=TRUE)
    mod8 <- try( fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
                           , hidden.areas=FALSE, trans.rate=trans.rate, assume.cladogenetic=TRUE) )
    saveRDS(mod8, file=paste0('ldg_geohisse', "_mod8.rds"))
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    capture.output( print( paste0("Analysis - ", 'ldg_geohisse', " - model 8 done.") ), file=paste0('ldg_geohisse',".log"), append=TRUE )
  }
  ## Model 3. Heterogeneous diversification, not tied to range evolution.
  ## Assumes three distinct diversification rates.
  if(model.number == 9){
    speciation <- c(1,1,1,2,2,2,3,3,3)
    extirpation <- c(1,1,2,2,3,3)
    trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=2, make.null=TRUE, include.jumps=FALSE, separate.extirpation=TRUE)
    mod9 <- try( fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
                           , hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE) )
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod9, file=paste0('ldg_geohisse', "_mod9.rds"))
  }
  ## Model 4. Heterogeneous diversification, tied to range evolution.
  ## Assumes 6 distinct diversification rates.
  if(model.number == 10){
    speciation <- c(1,2,3,4,5,6)
    extirpation <- c(1,2,3,4)
    trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=1, include.jumps=FALSE, separate.extirpation=TRUE)
    mod10 <- try( fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
                            , hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE) )
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod10, file=paste0('ldg_geohisse', "_mod10.rds"))
  }
  ## Model 5. Heterogeneous diversification, not tied to range evolution. Assumes 5 distinct diversification rates.
  if(model.number == 11){
    speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3))
    extirpation <- c(rep(1,2), rep(2,2), rep(3,2), rep(4,2), rep(5,2))
    trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=4, make.null=TRUE, include.jumps=FALSE, separate.extirpation=TRUE)
    mod11 <- try( fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
                            , hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE) )
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod11, file=paste0('ldg_geohisse', "_mod11.rds"))
  }
  ## Model 6. Heterogeneous diversification, not tied to range evolution. Assumes two distinct diversification rates.
  if(model.number == 12){
    speciation <- c(1,1,1,2,2,2)
    extirpation <- c(1,1,2,2)
    trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=1, make.null=TRUE, include.jumps=FALSE, separate.extirpation=TRUE)
    mod12 <- try( fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
                            , hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=TRUE) )
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod12, file=paste0('ldg_geohisse', "_mod12.rds"))
  }
  
  ###############################################################################
  ## Second block of anagenetic models MuSSE-like.
  ## These are very liberal models. They are really GeoSSE-like without cladogenetic events
  ###############################################################################
  
  ## Model 1. Transitions only. No character effect on diversification
  if(model.number == 13){
    speciation <- c(1,1,1)
    extirpation <- c(1,1,1)
    trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=0, include.jumps=FALSE, separate.extirpation=TRUE)
    mod13 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
                       , hidden.areas=FALSE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod13, file=paste0('ldg_geohisse', "_mod13.rds"))
    
  }
  ## Model 2. Character effect on diversification.
  if(model.number == 14){
    speciation <- c(1,2,3)
    extirpation <- c(1,2,3)
    trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=0, include.jumps=FALSE, separate.extirpation=TRUE)
    mod14 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
                       , hidden.areas=FALSE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod14, file=paste0('ldg_geohisse', "_mod14.rds"))
    
  }
  ## Model 3. No character effect on diversification.
  if(model.number == 15){
    speciation <- c(1,1,1,2,2,2,3,3,3)
    extirpation <- c(1,1,1,2,2,2,3,3,3)
    trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=2, include.jumps=FALSE, separate.extirpation=TRUE, make.null=TRUE)
    mod15 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
                       , hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod15, file=paste0('ldg_geohisse', "_mod15.rds"))
    
  }
  ## Model 4. Character effect on diversification, with a hidden state
  if(model.number == 16){
    speciation <- c(1,2,3,4,5,6)
    extirpation <- c(1,2,3,4,5,6)
    trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=1, include.jumps=FALSE, separate.extirpation=TRUE)
    mod16 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation, hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod16, file=paste0('ldg_geohisse', "_mod16.rds"))
    
  }
  ## Model 5. No character effect on diversification, multiple shifts
  if(model.number == 17){
    speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3))
    extirpation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3))
    trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=4, include.jumps=FALSE, separate.extirpation=TRUE, make.null=TRUE)
    mod17 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation, hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod17, file=paste0('ldg_geohisse', "_mod17.rds"))
    
  }
  ## Model 6*. No character effect on diversification, multiple shifts.
  if(model.number == 18){
    speciation <- c(rep(1,3), rep(2,3))
    extirpation <- c(rep(1,3), rep(2,3))
    trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=1, include.jumps=FALSE, separate.extirpation=TRUE, make.null=TRUE)
    mod18 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation, hidden.areas=TRUE
                       , trans.rate=trans.rate, assume.cladogenetic=FALSE)
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod18, file=paste0('ldg_geohisse', "_mod18.rds"))
    
  }
  ## Just trying to complete the sheeit -- four character independent shifts
  if(model.number == 19){
    speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3))
    extirpation <- c(rep(1,2), rep(2,2), rep(3,2), rep(4,2))
    trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=3, include.jumps=FALSE, separate.extirpation=FALSE, make.null=TRUE)
    mod19 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation, hidden.areas=TRUE
                       , trans.rate=trans.rate, assume.cladogenetic=TRUE)
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod19, file=paste0('ldg_geohisse', "_mod19.rds"))
    
  }
  if(model.number == 20){
    speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3))
    extirpation <- c(rep(1,2), rep(2,2), rep(3,2), rep(4,2))
    trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=3, include.jumps=FALSE, separate.extirpation=TRUE, make.null=TRUE)
    mod20 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation, hidden.areas=TRUE
                       , trans.rate=trans.rate, assume.cladogenetic=TRUE)
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod20, file=paste0('ldg_geohisse', "_mod20.rds"))
    
  }
  if(model.number == 21){
    speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3))
    extirpation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3))
    trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=3, include.jumps=FALSE, separate.extirpation=TRUE, make.null=TRUE)
    mod21 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation, hidden.areas=TRUE
                       , trans.rate=trans.rate, assume.cladogenetic=FALSE)
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod21, file=paste0('ldg_geohisse', "_mod21.rds"))
    
  }
  ## A special type of MuSSE model: No speciation or extinction in the widespread range. These are endemic only linked processes. We assume xA and xB are linked to transitions from AB.
  ## Model 1. Transitions only. No character effect on diversification
  if(model.number == 22){
    speciation <- c(1,1,0)
    extirpation <- c(1,1,0)
    trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=0, include.jumps=FALSE, separate.extirpation=FALSE)
    mod22 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
                       , hidden.areas=FALSE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod22, file=paste0('ldg_geohisse', "_mod22.rds"))
    
  }
  ## Model 2. Character effect on diversification.
  if(model.number == 23){
    speciation <- c(1,2,0)
    extirpation <- c(1,2,0)
    trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=0, include.jumps=FALSE, separate.extirpation=FALSE)
    mod23 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
                       , hidden.areas=FALSE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
    saveRDS(mod23, file=paste0('ldg_geohisse', "_mod23.rds"))
    
  }
  ## Model 3. No character effect on diversification.
  if(model.number == 24){
    speciation <- c(1,1,0,2,2,0,3,3,0)
    extirpation <- c(1,1,0,2,2,0,3,3,0)
    trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=2, include.jumps=FALSE, separate.extirpation=FALSE, make.null=TRUE)
    mod24 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
                       , hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod24, file=paste0('ldg_geohisse', "_mod24.rds"))
    
  }
  ## Model 4. Character effect on diversification, with a hidden state
  if(model.number == 25){
    speciation <- c(1,2,0,3,4,0)
    extirpation <- c(1,2,0,3,4,0)
    trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=1, include.jumps=FALSE, separate.extirpation=FALSE)
    mod25 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation, hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
    saveRDS(mod25, file=paste0('ldg_geohisse', "_mod25.rds"))
    
  }
  ## Model 5. No character effect on diversification, multiple shifts
  if(model.number == 26){
    speciation <- c(rep(1,2),0, rep(2,2),0, rep(3,2),0, rep(4,2),0, rep(5,2),0)
    extirpation <- c(rep(1,2),0, rep(2,2),0, rep(3,2),0, rep(4,2),0, rep(5,2),0)
    trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=4, include.jumps=FALSE, separate.extirpation=FALSE, make.null=TRUE)
    mod26 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation, hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod26, file=paste0('ldg_geohisse', "_mod26.rds"))
    
  }
  ## Model 6*. No character effect on diversification, multiple shifts.
  if(model.number == 27){
    speciation <- c(rep(1,2),0, rep(2,2),0)
    extirpation <- c(rep(1,2),0, rep(2,2),0)
    trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=1, include.jumps=FALSE, separate.extirpation=FALSE, make.null=TRUE)
    mod27 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation, hidden.areas=TRUE
                       , trans.rate=trans.rate, assume.cladogenetic=FALSE)
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod27, file=paste0('ldg_geohisse', "_mod27.rds"))
    
  }
  if(model.number == 28){
    speciation <- c(rep(1,2),0, rep(2,2),0, rep(3,2),0, rep(4,2),0)
    extirpation <- c(rep(1,2),0, rep(2,2),0, rep(3,2),0, rep(4,2),0)
    trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=3, include.jumps=FALSE, separate.extirpation=FALSE, make.null=TRUE)
    mod28 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation, hidden.areas=TRUE
                       , trans.rate=trans.rate, assume.cladogenetic=FALSE)
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod28, file=paste0('ldg_geohisse', "_mod28.rds"))
    
  }
  ## A special type of MuSSE model: No speciation or extinction in the widespread range. These are endemic only linked processes. We assume xA and xB are unlinked to transitions from AB.
  if(model.number == 29){
    speciation <- c(1,1,0)
    extirpation <- c(1,1,0)
    trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=0, include.jumps=FALSE, separate.extirpation=TRUE)
    mod29 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
                       , hidden.areas=FALSE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod29, file=paste0('ldg_geohisse', "_mod29.rds"))
    
  }
  ## Model 2. Character effect on diversification.
  if(model.number == 30){
    speciation <- c(1,2,0)
    extirpation <- c(1,2,0)
    trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=0, include.jumps=FALSE, separate.extirpation=TRUE)
    mod30 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
                       , hidden.areas=FALSE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod30, file=paste0('ldg_geohisse', "_mod30.rds"))
    
  }
  ## Model 3. No character effect on diversification.
  if(model.number == 31){
    speciation <- c(1,1,0,2,2,0,3,3,0)
    extirpation <- c(1,1,0,2,2,0,3,3,0)
    trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=2, include.jumps=FALSE, separate.extirpation=TRUE, make.null=TRUE)
    mod31 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation
                       , hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
    saveRDS(mod31, file=paste0('ldg_geohisse', "_mod31.rds"))
    
  }
  ## Model 4. Character effect on diversification, with a hidden state
  if(model.number == 32){
    speciation <- c(1,2,0,3,4,0)
    extirpation <- c(1,2,0,3,4,0)
    trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=1, include.jumps=FALSE, separate.extirpation=TRUE)
    mod32 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation, hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod32, file=paste0('ldg_geohisse', "_mod32.rds"))
    
  }
  ## Model 5. No character effect on diversification, multiple shifts
  if(model.number == 33){
    speciation <- c(rep(1,2),0, rep(2,2),0, rep(3,2),0, rep(4,2),0, rep(5,2),0)
    extirpation <- c(rep(1,2),0, rep(2,2),0, rep(3,2),0, rep(4,2),0, rep(5,2),0)
    trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=4, include.jumps=FALSE, separate.extirpation=TRUE, make.null=TRUE)
    mod33 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation, hidden.areas=TRUE, trans.rate=trans.rate, assume.cladogenetic=FALSE)
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod33, file=paste0('ldg_geohisse', "_mod33.rds"))
    
  }
  ## Model 6*. No character effect on diversification, multiple shifts.
  if(model.number == 34){
    speciation <- c(rep(1,2),0, rep(2,2),0)
    extirpation <- c(rep(1,2),0, rep(2,2),0)
    trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=1, include.jumps=FALSE, separate.extirpation=TRUE, make.null=TRUE)
    mod34 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation, hidden.areas=TRUE
                       , trans.rate=trans.rate, assume.cladogenetic=FALSE)
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod34, file=paste0('ldg_geohisse', "_mod34.rds"))
    
  }
  if(model.number == 35){
    speciation <- c(rep(1,2),0, rep(2,2),0, rep(3,2),0, rep(4,2),0)
    extirpation <- c(rep(1,2),0, rep(2,2),0, rep(3,2),0, rep(4,2),0)
    trans.rate <- trans.rate <- TransMatMakerfGeoHiSSE(hidden.areas=3, include.jumps=FALSE, separate.extirpation=TRUE, make.null=TRUE)
    mod35 <- fGeoHiSSE(phy, dat, f=c(sampling0,sampling1,sampling2), turnover=speciation, extinct.frac=extirpation, hidden.areas=TRUE
                       , trans.rate=trans.rate, assume.cladogenetic=FALSE)
    capture.output( print( paste0("Analysis - ", model.number, "done.") ), file=paste0("ldg_geohisse_",model.number,".log") )
    saveRDS(mod35, file=paste0('ldg_geohisse', "_mod35.rds"))
    
  }
  
  
}

## Load data and evaluate models. Here we will make 80 replicates for dataset. This will maximize the use of the cores over time.
args<-commandArgs(trailingOnly = TRUE)
print(args)

geohisse.ldgplants.evaluate.models(model.number=args[1])





