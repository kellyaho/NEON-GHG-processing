#####################################################################
# Rheadspace_GHG.R
#
# R function that calculates dissolved CO2, CH4, and N2O partial pressures and 
# concentration from pre- and post-equilibration headspace mixing ratios.  This 
# function is modeled after Rheadspace.R by Marce, Kim, and Prairie. 
#
# This function calculates three versions of dissolved GHG outputs using three 
# different pre-equilbiration headspace mixing ratios : 1. mixing ratio of paired 
# ATM sample, 2. estimate calculated from loess smoothing of the ATM samples by site, 
# 3. median atmospheric mixing value by site. These three pre-equilibration headspaces
# can be calculated using the create_input.R script.
#
# For CO2, the calculation accounts for carbonate equilbrium according to 
# Rheadspace.R by Marce, Kim, and Prairie as in 
# Koschorreck, M., Y.T. Prairie, J. Kim, and R. Marcé. 2020. 
# Technical note: CO2 is not like CH4 – limits of the headspace method to analyse pCO2 in water. Biogeosciences, 2021
# and https://github.com/icra/headspace/blob/master/Rheadspace.R.
# This function modifies the Rheadspace function by using only freshwater constants  
# and updates the Henry’s law constants for CO2 according to Sander 2015.
#
# This function also calculates dissolved CH4 and N2O using CH4 and N2O Henry's 
# Law constants according to Sander 2015.
#######################################################################
#######
Rheadspace_GHG <-  function(...){
  arguments <- list(...)
  
  if (is.data.frame(arguments[[1]])) {
    input.table=arguments[[1]]
    if (dim(input.table)[2]!=21){
      stop("Input dataframe error", call.=FALSE)
    }else{
      Sample.ID = as.character(input.table$Sample.ID)
      timestamp = input.table$timestamp
      mCH4_headspace.paired = input.table$HS.mCH4.before #the mCH4 (ppmv) of the headspace "before" equilibration - paired
      mCO2_headspace.paired = input.table$HS.mCO2.before #the mCO2 (ppmv) of the headspace "before" equilibration - paired
      mN2O_headspace.paired = input.table$HS.mN2O.before #the mN2O (ppmv) of the headspace "before" equilibration - paired
      mCH4_headspace.smoothed = input.table$smoothed.HS.mCH4.before #the mCH4 (ppmv) of the headspace "before" equilibration - smoothed
      mCO2_headspace.smoothed = input.table$smoothed.HS.mCO2.before #the mCO2 (ppmv) of the headspace "before" equilibration - smoothed
      mN2O_headspace.smoothed = input.table$smoothed.HS.mN2O.before #the mN2O (ppmv) of the headspace "before" equilibration - smoothed
      mCH4_headspace.median = input.table$median.HS.mCH4.before #the mCH4 (ppmv) of the headspace "before" equilibration - site median
      mCO2_headspace.median = input.table$median.HS.mCO2.before #the mCO2 (ppmv) of the headspace "before" equilibration - site median
      mN2O_headspace.median = input.table$median.HS.mN2O.before #the mN2O (ppmv) of the headspace "before" equilibration - site median
      mCH4_eq = input.table$HS.mCH4.after #the measured mCH4 (ppmv) of the headspace "after" equilibration
      mCO2_eq = input.table$HS.mCO2.after #the measured mCO2 (ppmv) of the headspace "after" equilibration
      mN2O_eq = input.table$HS.mN2O.after #the measured mN2O (ppmv) of the headspace "after" equilibration
      temp_insitu = input.table$Temp.insitu #in situ water temperature in degrees celsius
      temp_eq = input.table$Temp.equil #the water temperature after equilibration in degree celsius
      alk = input.table$Alkalinity.measured #Total alkalinity (micro eq/L) of the water sample
      vol_gas = input.table$Volume.gas #Volume of gas in the headspace vessel (mL)
      vol_water = input.table$Volume.water #Volume of water in the headspace vessel (mL)   
      Bar.pressure = input.table$Bar.pressure #Barometric pressure at field conditions in kPa. 101.325 kPa = 1 atm   
      site = input.table$site
    } 
  } else if (length(arguments)==13) {
    Sample.ID = as.character(arguments[[1]])
    timestamp = arguments[[2]]
    mCH4_headspace.paired = arguments[[3]] #the mCH4 (ppmv) of the headspace "before" equilibration - paired
    mCO2_headspace.paired = arguments[[4]] #the mCO2 (ppmv) of the headspace "before" equilibration - paired
    mN2O_headspace.paired = arguments[[5]] #the mN2O (ppmv) of the headspace "before" equilibration - paired
    mCH4_headspace.smoothed = arguments[[6]] #the mCH4 (ppmv) of the headspace "before" equilibration - smoothed
    mCO2_headspace.smoothed = arguments[[7]] #the mCO2 (ppmv) of the headspace "before" equilibration - smoothed
    mN2O_headspace.smoothed = arguments[[8]] #the mN2O (ppmv) of the headspace "before" equilibration - smoothed
    mCH4_headspace.median = arguments[[9]] #the mCH4 (ppmv) of the headspace "before" equilibration - median
    mCO2_headspace.median = arguments[[10]] #the mCO2 (ppmv) of the headspace "before" equilibration - median
    mN2O_headspace.median = arguments[[11]] #the mN2O (ppmv) of the headspace "before" equilibration - median
    mCH4_eq = arguments[[12]] #the measured mCH4 (ppmv) of the headspace "after" equilibration
    mCO2_eq = arguments[[13]] #the measured mCO2 (ppmv) of the headspace "after" equilibration
    mN2O_eq = arguments[[14]] #the measured mN2O (ppmv) of the headspace "after" equilibration
    temp_insitu = arguments[[15]] #in situ water temperature in degrees celsius
    temp_eq = arguments[[16]] #the water temperature after equilibration in degree celsius
    alk = arguments[[17]] #Total alkalinity (micro eq/L) of the water sample
    vol_gas = arguments[[18]] #Volume of gas in the headspace vessel (mL)
    vol_water = arguments[[19]] #Volume of water in the headspace vessel (mL)   
    Bar.pressure = arguments[[20]] #Barometric pressure at field conditions in kPa. 101.325 kPa = 1 atm   
    site = arguments[[21]]
  } else {
    stop("Input error", call.=FALSE)
  }
  
  pGHG_orig <- data.frame(matrix(NA,length(mCO2_headspace.paired),19))
  names(pGHG_orig) <- c("Sample.ID",
                        "pCH4.uatm.paired", 
                        "CH4.uM.paired",
                        "pCH4.uatm.smoothed", 
                        "CH4.uM.smoothed",
                        "pCH4.uatm.median", 
                        "CH4.uM.median",
                        "pCO2.uatm.paired", 
                        "CO2.uM.paired", 
                        "pCO2.uatm.smoothed", 
                        "CO2.uM.smoothed", 
                        "pCO2.uatm.median", 
                        "CO2.uM.median", 
                        "pN2O.uatm.paired", 
                        "N2O.uM.paired",
                        "pN2O.uatm.smoothed", 
                        "N2O.uM.smoothed",
                        "pN2O.uatm.median", 
                        "N2O.uM.median")
  
  R <- 0.082057338 #L atm K-1 mol-1
  
  for (i in 1:length(mCO2_headspace.paired)){ 
    
    AT = alk[i]*(1e-6) #ueq/l to mol/L
    
    #Constants of the carbonate equilibrium
    # Kw = the dissociation constant of H2O into H+ and OH- (Dickson and Riley, 1979)
    # K1 = the equilibrium constant between CO2 and HCO3- (Millero 1979)
    # K2 = the equilibrium constant between HCO3- and CO3 2- (Millero 1979)
    
    # Kh.co2 = the solubility of CO2 in water - equilibration conditions
    # Kh2.co2 = the solubility of CO2 in water - in situ field conditions
    # Kh.ch4 = the solubility of CH4 in water - equilibration conditions
    # Kh2.ch4 = the solubility of CH4 in water - in situ field conditions
    # Kh.n2o = the solubility of N2O in water - equilibration conditions
    # Kh2.n2o = the solubility of N2O in water - in situ field conditions
    # Solubility coefficients from Sander 2015 (1.4 x 10^-5 mol m^-3 Pa ^-1 and 
    # 1900 K for CH4, 3.3 x 10^-4 mol m^-3 Pa ^-1 and 2400 K for CO2, and 
    # 2.4 x 10^-4 mol m^-3 Pa ^-1 and 2700 K for N2O)

    K1=10^-(-126.34048+6320.813/(temp_eq[i]+273.15)+19.568224*log(temp_eq[i]+273.15))
    K2=10^-(-90.18333+5143.692/(temp_eq[i]+273.15)+14.613358*log(temp_eq[i]+273.15))
    Kw = exp(148.9652-13847.26/(temp_eq[i]+273.15)-23.6521*log(273.15+temp_eq[i]))
    
    Kh.co2 = 0.00033*exp(2400*(1/(273.15+temp_eq[i])-1/298.15))*101325/1000 # mol/L/atm equilibration conditions
    Kh2.co2 = 0.00033*exp(2400*(1/(273.15+temp_insitu[i])-1/298.15))*101325/1000 # mol/L/atm original conditions
    Kh.ch4 = 0.000014*exp(1900*(1/(273.15+temp_eq[i])-1/298.15))*101325/1000# mol/L/atm equilibration conditions
    Kh2.ch4 = 0.000014*exp(1900*(1/(273.15+temp_insitu[i])-1/298.15))*101325/1000 # mol/L/atm original conditions
    Kh.n2o = 0.00024*exp(2700*(1/(273.15+temp_eq[i])-1/298.15))*101325/1000 # mol/L/atm equilibration conditions
    Kh2.n2o = 0.00024*exp(2700*(1/(273.15+temp_insitu[i])-1/298.15))*101325/1000 # mol/L/atm original conditions
    
    HS.ratio <- vol_gas[i]/vol_water[i] 
    
    ##### CH4 ####
    CH4_solution <- mCH4_eq[i]/1000000*Kh.ch4 #mol/L
    CH4_solution_mass <- CH4_solution * vol_water[i]/1000 #mol
    
    final_C_headspace_mass.ch4 <- mCH4_eq[i]/1000000*(vol_gas[i]/1000) / (R * (temp_eq[i]+273.15)) #mol
    
    mols_headspace.ch4.paired <- mCH4_headspace.paired[i]/1000000*(vol_gas[i]/1000)/(R * (temp_eq[i]+273.15)) #mol PV / RT = n
    mols_headspace.ch4.smoothed <- mCH4_headspace.smoothed[i]/1000000*(vol_gas[i]/1000)/(R * (temp_eq[i]+273.15)) #mol PV / RT = n
    mols_headspace.ch4.median <- mCH4_headspace.median[i]/1000000*(vol_gas[i]/1000)/(R * (temp_eq[i]+273.15)) #mol PV / RT = n
    Sample_CH4_mass.paired <- CH4_solution_mass + final_C_headspace_mass.ch4 - mols_headspace.ch4.paired #mol
    Sample_CH4_mass.smoothed <- CH4_solution_mass + final_C_headspace_mass.ch4 - mols_headspace.ch4.smoothed #mol
    Sample_CH4_mass.median <- CH4_solution_mass + final_C_headspace_mass.ch4 - mols_headspace.ch4.median #mol
    Sample_CH4_conc.paired <- Sample_CH4_mass.paired/(vol_water[i]/1000) #mol/L
    Sample_CH4_conc.smoothed <- Sample_CH4_mass.smoothed/(vol_water[i]/1000) #mol/L
    Sample_CH4_conc.median <- Sample_CH4_mass.median/(vol_water[i]/1000) #mol/L
    
    pGHG_orig[i,1] <- as.character(Sample.ID[i])
    pGHG_orig[i,2] <- Sample_CH4_conc.paired/Kh2.ch4*1000000*Bar.pressure[i]/101.325 # micro-atm
    pGHG_orig[i,3] <- Sample_CH4_conc.paired*1000000 # micro-mol/L
    pGHG_orig[i,4] <- Sample_CH4_conc.smoothed/Kh2.ch4*1000000*Bar.pressure[i]/101.325 # micro-atm
    pGHG_orig[i,5] <- Sample_CH4_conc.smoothed*1000000 # micro-mol/L
    pGHG_orig[i,6] <- Sample_CH4_conc.median/Kh2.ch4*1000000*Bar.pressure[i]/101.325 # micro-atm
    pGHG_orig[i,7] <- Sample_CH4_conc.median*1000000 # micro-mol/L
    
    ###CO2#### Adapted from Rheadspace.R by Marce, Kim, and Prairie
    #DIC at equilibrium
    co2 <- Kh.co2 * mCO2_eq[i]/1000000
    h_all <- polyroot(c(-(2*K1*K2*co2),-(co2*K1+Kw),AT,1))
    real<-Re(h_all)
    h <-real[which(real>0)]
    DIC_eq <- co2 * (1 + K1/h + K1 * K2/(h * h))
    
    #DIC in the original sample
    DIC_ori.paired <- DIC_eq + (mCO2_eq[i] - mCO2_headspace.paired[i])/1000000/(R*(temp_eq[i]+273.15))*HS.ratio
    DIC_ori.smoothed <- DIC_eq + (mCO2_eq[i] - mCO2_headspace.smoothed[i])/1000000/(R*(temp_eq[i]+273.15))*HS.ratio
    DIC_ori.median <- DIC_eq + (mCO2_eq[i] - mCO2_headspace.median[i])/1000000/(R*(temp_eq[i]+273.15))*HS.ratio
    
    #pCO2 in the original sample
    h_all.paired <- polyroot(c(-(K1*K2*Kw),K1*K2*AT-K1*Kw-2*DIC_ori.paired*K1*K2,AT*K1-Kw+K1*K2-DIC_ori.paired*K1,AT+K1,1))
    h_all.smoothed <- polyroot(c(-(K1*K2*Kw),K1*K2*AT-K1*Kw-2*DIC_ori.smoothed*K1*K2,AT*K1-Kw+K1*K2-DIC_ori.smoothed*K1,AT+K1,1))
    h_all.median <- polyroot(c(-(K1*K2*Kw),K1*K2*AT-K1*Kw-2*DIC_ori.median*K1*K2,AT*K1-Kw+K1*K2-DIC_ori.median*K1,AT+K1,1))
    
    real.paired<-Re(h_all.paired)
    real.smoothed<-Re(h_all.smoothed)
    real.median<-Re(h_all.median)
    h.paired <-real.paired[which(real.paired>0)]
    h.smoothed <-real.smoothed[which(real.smoothed>0)]
    h.median <-real.median[which(real.median>0)]
    
    co2.paired <- h.paired* (DIC_ori.paired * h.paired * K1/(h.paired * h.paired + K1 * h.paired + K1 * K2)) / K1
    co2.smoothed <- h.smoothed* (DIC_ori.smoothed * h.smoothed * K1/(h.smoothed * h.smoothed + K1 * h.smoothed + K1 * K2)) / K1
    co2.median <- h.median* (DIC_ori.median * h.median * K1/(h.median * h.median + K1 * h.median + K1 * K2)) / K1
    
    pGHG_orig[i,8] <- co2.paired/Kh2.co2*1000000*Bar.pressure[i]/101.325
    pGHG_orig[i,9] <- co2.paired*1000000
    pGHG_orig[i,10] <- co2.smoothed/Kh2.co2*1000000*Bar.pressure[i]/101.325
    pGHG_orig[i,11] <- co2.smoothed*1000000
    pGHG_orig[i,12] <- co2.median/Kh2.co2*1000000*Bar.pressure[i]/101.325
    pGHG_orig[i,13] <- co2.median*1000000
    
    #N2O####
    N2O_solution <- mN2O_eq[i]/1000000*Kh.n2o #mol/L
    N2O_solution_mass <- N2O_solution * vol_water[i]/1000 #mol
    
    final_N_headspace_mass <- mN2O_eq[i]/1000000*(vol_gas[i]/1000) / (R * (temp_eq[i]+273.15)) #mol
    
    mols_headspace.n2o.paired <- mN2O_headspace.paired[i]/1000000*(vol_gas[i]/1000)/(R * (temp_eq[i]+273.15)) #mol PV / RT = n
    mols_headspace.n2o.smoothed <- mN2O_headspace.smoothed[i]/1000000*(vol_gas[i]/1000)/(R * (temp_eq[i]+273.15)) #mol PV / RT = n
    mols_headspace.n2o.median <- mN2O_headspace.median[i]/1000000*(vol_gas[i]/1000)/(R * (temp_eq[i]+273.15)) #mol PV / RT = n
    
    Sample_N2O_mass.paired <- N2O_solution_mass + final_N_headspace_mass - mols_headspace.n2o.paired #mol
    Sample_N2O_mass.smoothed <- N2O_solution_mass + final_N_headspace_mass - mols_headspace.n2o.smoothed #mol
    Sample_N2O_mass.median <- N2O_solution_mass + final_N_headspace_mass - mols_headspace.n2o.median #mol
    Sample_N2O_conc.paired <- Sample_N2O_mass.paired/(vol_water[i]/1000) #mol/L
    Sample_N2O_conc.smoothed <- Sample_N2O_mass.smoothed/(vol_water[i]/1000) #mol/L
    Sample_N2O_conc.median <- Sample_N2O_mass.median/(vol_water[i]/1000) #mol/L
    
    pGHG_orig[i,14] <- Sample_N2O_conc.paired/Kh2.n2o*1000000*Bar.pressure[i]/101.325 # micro-atm
    pGHG_orig[i,15] <- Sample_N2O_conc.paired*1000000 # micro-mol/L
    pGHG_orig[i,16] <- Sample_N2O_conc.smoothed/Kh2.n2o*1000000*Bar.pressure[i]/101.325 # micro-atm
    pGHG_orig[i,17] <- Sample_N2O_conc.smoothed*1000000 # micro-mol/L
    pGHG_orig[i,18] <- Sample_N2O_conc.median/Kh2.n2o*1000000*Bar.pressure[i]/101.325 # micro-atm
    pGHG_orig[i,19] <- Sample_N2O_conc.median*1000000 # micro-mol/L
  }
  
  
  return(pGHG_orig) 
  
}
