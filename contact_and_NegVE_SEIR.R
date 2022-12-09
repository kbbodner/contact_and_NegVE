###################################################################### #
#                                                                       #                                                   #
#  Vaccinated Contact Heterogeneity and Negative Vaccine Effectiveness  #
#                       Main File (SEIR Model Results)                  #
#                                                                       #
########################################################################

#PACKAGES
require(deSolve)
require(tidyverse)
require(ggplot2)
require(reshape)
library(ggpubr)
library(plotly)


###DATA GENERATION FOR FIGURE 1 AND SUPPLEMENTARY MATERIAL WEB FIGURE 3 

#VE function
#TWO INPUTS:
  #vaccine.type is either "AllOrNothing" or "Leaky", default is "AllOrNothing"
  #contact.type is either "Homog." or "Heterog.", default is "Homog." 
      #(Heterog. is vaccinated contact heterogeneity with 50% higher contact among vaccinated)
#OUTPUT: Dataset of VE measurements and prop infected from an SEIR model using 175 timesteps
#Dataset contains simulations for VE_S (vaccine efficacy against susc.) and VE_I (vaccine efficacy against infect.) at 0.1, 0.5

VE_function<-function(vaccine.type="AllOrNothing", contact.type="Homog."){
  
  
  names<-c("time","S1","S2","E1","E2","I1","I2","R1","R2",
           "CumIncid1","CumIncid2","Incid1","Incid2",
           "a","b","c", "d","OR.overall",
           "RR.overall","prop.unvac.overall","prop.vac.overall", "prop.total.overall",
           "e","f","g","h", "OR.perday","RR.perday","prop.unvac.perday","prop.vac.perday",
           "label.ve.susc","label.ve.infect")
  
  #Empty df for all VE combos
  Main.VE.df<-data.frame(matrix(ncol=length(names),nrow=0))
  
  #Add column names to df
  colnames(Main.VE.df)<-names
  
  #Full ODE
  ve_model<-function(t,x,parms){ 
    
    
    ncompartment <- 6 # 6 compartments (S, E, I, R)
    ngroup <- length(x)/ncompartment   #number of groups
    S      <- as.matrix(x[1:ngroup])                          #compartment 1
    E    <- as.matrix(x[(ngroup+1):(2*ngroup)])               #compartment 2
    I    <- as.matrix(x[(2*ngroup+1):(3*ngroup)])             #compartment 3
    R   <- as.matrix(x[(3*ngroup+1):(4*ngroup)])              #compartment 4
    CumIncid    <- as.matrix(x[(4*ngroup+1):(5*ngroup)])      #compartment 5
    Incid   <- as.matrix(x[(5*ngroup+1):(6*ngroup)])          #compartment 6
    
    with(as.list(parms),{
      
      N <- S+E+I+R #this just creates a proportional metric per group    
      rho=c(beta,beta*alpha) 
      sigma.matrix=matrix(c(1,1,sigma,sigma),nrow=2,ncol=2)
      
      
      # ODEs
      
      dS <- -as.matrix(S)*rho*(as.matrix(C)*sigma.matrix)%*%as.matrix(I/N)
      
      dE <- as.matrix(S)*rho*(as.matrix(C)*sigma.matrix)%*%as.matrix(I/N) - mu*as.matrix(E) #we assume exposed state isn't infectious
      
      dI <- +mu*as.matrix(E) - gamma*as.matrix(I)       
      
      dR <- +gamma*as.matrix(I)
      
      dCumIncid <-+as.matrix(S)*rho*(as.matrix(C)*sigma.matrix)%*%as.matrix(I/N)
      
      dIncid <-+as.matrix(S)*rho*(as.matrix(C)*sigma.matrix)%*%as.matrix(I/N)  - as.matrix(Incid)
      
      
      #the output has to be in the same order as the model compartments in the begining of the function
      out = c(dS,
              dE,
              dI,         
              dR,          
              dCumIncid,
              dIncid)
      list(out)
      
    }) 
  }
  
  ###General Model Set-up (Applicable to all Model Runs)
  #SET-UP MODEL
  # initial conditions
  npop   <- 100000            # 100,000 total population 
  f      <- c(0.25, 0.75)     # 2 subgroups f[1]=  unvaccinated; f[2] = vaccinated
  Ntot   <- npop*f            # population size in each subgroup (vector of 2 elements)
  ngroup <- length(f)         # length of vector representing number of subgroups
  I_init    <- c(1,1)         # 1 infectious individual per group
  E_init    <- c(0,0)         # 0 exposed individuals
  
  
  cr  = c(6,6)       #baseline number of contacts per day among unvac[1] and vack[2]
  epsilon = 1      #where 1 = proportionate; and 0 = assortative (who contacts whom)
  
  rho_unvac_unvac = (1-epsilon) + epsilon*(Ntot[1]*cr[1]/(Ntot[1]*cr[1] + Ntot[2]*cr[2]))
  rho_unvac_vac  =               epsilon*(Ntot[2]*cr[2]/(Ntot[1]*cr[1] + Ntot[2]*cr[2]))
  rho_vac_unvac   =             epsilon*(Ntot[1]*cr[1]/(Ntot[1]*cr[1] + Ntot[2]*cr[2]))
  rho_vac_vac  =(1-epsilon) + epsilon*(Ntot[2]*cr[2]/(Ntot[1]*cr[1] + Ntot[2]*cr[2]))
  
  weighted_average_cr = (cr[1]*Ntot[1] + cr[2]*Ntot[2])/(Ntot[1]+Ntot[2])
  
  
  #Homogenius or Heterogeneous
  C   <- matrix(0,nrow=ngroup,ncol=ngroup) # contact/mixing matrix of ngroup by ngroup
  
  contact.increase<-1.5 #only bumps you up to 8.25 compared to 6
  
  if (contact.type=="Homog."){
    C[1,1]  <- cr[1]*rho_unvac_unvac  # contacts between unvac and unvac
    C[1,2]  <- cr[1]*rho_unvac_vac   # contacts between unvac and vac
    C[2,1]  <- cr[2]*rho_vac_unvac   # contacts between vac and unvac
    C[2,2]  <- cr[2]*rho_vac_vac    # contacts between vac and vac
  }else if (contact.type=="Heterog."){
    C[1,1]  <- cr[1]*rho_unvac_unvac # contacts between unvac and unvac
    C[1,2]  <- cr[1]*rho_unvac_vac   # contacts between unvac and vac
    C[2,1]  <- cr[2]*rho_vac_unvac   # contacts between vac and unvac
    C[2,2]  <- cr[2]*rho_vac_vac*contact.increase   # contacts between vac and vac
  }
  
  #specify other parameters
  gamma   <- 1/10 #recovery rate
  mu<- 1/4 #latency period #2 to 4 days (or 3-4)
  beta<-0.1 #this is for an R0 of 6 
  
  
  #VACCINE EFFICACY ITERATIONS
  ve.start<-0.1
  ve.susc<-ve.start
  
  #Vaccine Efficacy Against Infectiousness Loop
  for (i in 1:2){
    
    
    ve.infect<-ve.start
    
    #Vaccine Efficacy Against Susceptibility Loop
    for(j in 1:2){
      
      
      sigma<-1-ve.infect 
      
      #Here we incorporate VE against suscept. depending on Vaccine Type (all or nothing or leaky)
      if(vaccine.type=="AllOrNothing"){
        alpha=(1-0) 
        R_init<- c(0,Ntot[2]*ve.susc)  
        
      }else if(vaccine.type=="Leaky"){
        alpha=(1-ve.susc)
        R_init<- c(0,0)  
      }
      
      #ODE INitial Conditions & Pars
      S_init    <- Ntot - E_init - I_init - R_init   # initial number susceptible, subtracting number of exposed, infectious, and immune from in each group
     
      rho=c(beta,beta*alpha) 
      sigma.matrix=matrix(c(1,1,sigma,sigma),nrow=2,ncol=2)
      
      CumIncid_init <- I_init
      Incid_init <- I_init
      
      parm <-c(gamma=gamma,
               mu=mu,
               beta=beta, 
               alpha=alpha,
               sigma=sigma,
               C=C)
      
      init <-c(S = S_init,
               E = E_init,
               I = I_init,           
               R = R_init,               
               CumIncid = CumIncid_init,
               Incid = Incid_init)  
      
      
      time = seq(0,175,1)  #run model for 175 days
      
      #RUN MODEL
      
      model.output<-as.data.frame(lsoda(y=init,times=time,func=ve_model,parms=parm)) 
      
      
      
      #Calculating Cumulative Incidence OR and RR
      model.output$a<-model.output$CumIncid2
      model.output$b<-Ntot[2] -model.output$CumIncid2 
      model.output$c<-model.output$CumIncid1
      model.output$d<-Ntot[1] -model.output$CumIncid1
      
      model.output$OR.overall<-1-((model.output$a*model.output$d)/(model.output$b*model.output$c))
      #model.output$RR.overall<-1-((model.output$a/(model.output$a + model.output$b))/(model.output$c/(model.output$c + model.output$d)))
      model.output$RR.overall<-1-((model.output$a/Ntot[2])/(model.output$c/Ntot[1]))
      
      #Calculating Cumulative Incidence Proportion Infected
      model.output$prop.unvac.overall<-model.output$CumIncid1/Ntot[1]
      model.output$prop.vac.overall<-model.output$CumIncid2/Ntot[2]
      model.output$prop.total.overall<-(model.output$CumIncid1 + model.output$CumIncid2)/sum(Ntot)
      
      #Calculating Incidence OR and RR
      model.output$e<-model.output$Incid2
      model.output$f<-Ntot[2] -model.output$Incid2
      model.output$g<-model.output$Incid1
      model.output$h<-Ntot[1] -model.output$Incid1
      
      model.output$OR.perday<-1 - ((model.output$e*model.output$h)/(model.output$f*model.output$g))
      model.output$RR.perday<-1-((model.output$e/(model.output$e + model.output$f))/(model.output$g/(model.output$g + model.output$h)))
      
      #Calculating Incidence Proportion Infected
      model.output$prop.unvac.perday<-model.output$I1/Ntot[1]
      model.output$prop.vac.perday<-model.output$I2/Ntot[2]
      
      #create variable indicator for VEs
      model.output$label.ve.susc=paste("VE Susc.=", ve.susc,sep=" ")
      model.output$label.ve.infect=paste("VE Infect.=", ve.infect, sep=" ")
      
      #Add results to main dataset
      Main.VE.df<-merge(Main.VE.df,model.output,all.x=T,all.y=T,sort=F)
      
      ve.infect<-ve.infect + 0.4
    }
    ve.susc<-ve.susc + 0.4
    
  }
  
  #Add to Vaccine Type and whether contact is Heterog or Homog.
  Main.VE.df$vaccine.type<-paste(vaccine.type)
  Main.VE.df$contact.type<-paste(contact.type)
  
  
  return(Main.VE.df)
  
}


#Create datasets using VE_function
AoN_Homog_VE<-VE_function(vaccine.type="AllOrNothing", contact.type="Homog.")
AoN_Heterog_VE<-VE_function(vaccine.type="AllOrNothing", contact.type="Heterog.")

#Final dataset
subset.homhet<-rbind(AoN_Homog_VE, AoN_Heterog_VE)

#Create label grouping with VE against Susc
subset.homhet$indicator.contact.susc<-paste(subset.homhet$contact.type, subset.homhet$label.ve.susc, sep="; ")

#Create overall proportion infected (exposed and infectious are both considered infected)
subset.homhet$prop.infected<-(subset.homhet$E1 + subset.homhet$E2 +subset.homhet$I1 + subset.homhet$I2)/(subset.homhet$E1 + subset.homhet$E2 + subset.homhet$I1 + subset.homhet$I2 + subset.homhet$S1 + subset.homhet$S2 + subset.homhet$R1 + subset.homhet$R2)


###FIGURE 1 (a and b) LINE PLOTS (VE based on Relative Risk and Proportion Infected)

colours.hom.het<-c("#56B4E9", "#0072B2", "#E69F00", "#D55E00")


#Vaccine Effectiveness (1-RR) Plot
#Note that RR.overall is vaccine effectiveness calculated using the relative risk
#t=51 and t=75 when RR crosses over in heterog case with low VE to be positive (see Supplementary Material Web Appendix 3 section below)
Ve.RR.overall.plot<-ggplot(subset.homhet, aes(x=time, y=RR.overall, color=indicator.contact.susc, linetype=label.ve.infect)) +
  theme_classic() + geom_hline(yintercept=0, linetype="dotted", color = "darkgrey", size=0.8) +
  geom_segment(aes(x=51, y=-Inf, xend=51, yend=0),linetype="dotted", color = "darkgrey", size=0.8) +
  geom_segment(aes(x=75, y=-Inf, xend=75, yend=0),linetype="dotted", color = "darkgrey", size=0.8) +
  geom_line(size=1.5) +
  ggtitle("Vaccine Effectiveness and Proportion Infected Across Time") + 
  labs(y="Vaccine Effectiveness (1-Relative Risk)", x = "Days", color="VE Values") + 
  scale_x_continuous(limits = c(0, 175),breaks=seq(0,175,25) ) + scale_y_continuous(limits = c(-0.3, 0.8)) +
  scale_color_manual(values=colours.hom.het) + theme(axis.title.x=element_blank(),
                                                     axis.text.x=element_blank(),
                                                     axis.ticks.x=element_blank())

#Prop.Infected
prop.infected.plot<-ggplot(subset.homhet, aes(x=time, y=prop.infected, color=indicator.contact.susc, linetype=label.ve.infect)) +
  theme_classic() +
  geom_line(size=1.5) +
  labs(y="Proportion Infected", x = "Days", color="VE Values") + 
  scale_x_continuous(limits = c(0, 175), breaks=seq(0,175,25) ) +
  scale_color_manual(values=colours.hom.het) 


line.plots.overall<-ggarrange(Ve.RR.overall.plot, prop.infected.plot, ncol=1, common.legend=T, legend = "right")

line.plots.overall #Figure 1


###FIGURE 2 CONTOUR PLOTS 
#Figure 2a: Vaccine Efficacy against Susceptibility (VES) with Vaccine Efficacy against Infectiousness (VEI) 
#Figure 2b: Vaccine Efficacy against Susceptibility (VES) with level of vaccinated contact heterogeneity

#Read-in data for contour plots (data generated from contour_data_generation_SEIR.R)
VeT.contour.df<-read.csv("data/VES_and_VEI_SEIR.csv" )
contact.VES.df<-read.csv("data/VES_and_ContactIncrease_SEIR.csv")

#Subsetting VES<0.51 as higher values produces the same results as 0.5
veT.susc0.5<-subset(VeT.contour.df, ve.susc<0.51)

#Contour Plot for VES and VEI
fig.veT<-plot_ly(z = veT.susc0.5$min.RR.overall,x=veT.susc0.5$ve.susc, 
                 y=veT.susc0.5$ve.infect, type = "contour", autocontour=F,line = list(width = 1, color = "black"),
                 contours = list(coloring = 'heatmap',  start = -.65, end = 0, size = 0.05, showlabels=T, 
                                 labelfont = list(size = 27)),
                 colorscale='Jet',reversescale=T,
                 line = list(smoothing = 0.3))
fig.veT <- fig.veT %>% colorbar(title = "Vaccine \nEffectiveness (1-RR)", titlefont = list(size = 15))  
fig.veT <- fig.veT%>% layout(xaxis = list(title = "VE against susceptibility"), 
                             yaxis = list(title = "VE against infectiousness")) 


fig.veT

#Save as .eps (orca must be installed) [it says deprecated but still is operational]
#orca(fig.veT, "VET_contour", format="eps", width=1200)

#Subsetting VES<0.51 as higher values produces the same results as 0.5
contact.ve.susc0.5<-subset(contact.VES.df, ve.susc<0.51)

#Convert contact rate into % increase above Homogeneous contact scenario (i.e. random mixing scenario)
contact.ve.susc0.5$contact.increase.percent<-(contact.ve.susc0.5$contact.increase - 1)*100

#Contour Plot for VES and contact heterogeneity
fig.contact <-plot_ly(z = contact.ve.susc0.5$min.RR.overall,x=contact.ve.susc0.5$ve.susc, 
                      y=contact.ve.susc0.5$contact.increase.percent,
                      type = "contour", autocontour=F, line = list(width = 1, color = "black"),
                      contours = list(coloring = 'heatmap',  start = -.65, end = 0, size = 0.05, showlabels=T,
                                      labelfont = list(size = 27)),
                      colorscale='Jet',reversescale=T,
                      line = list(smoothing = 0.3))
#fig.contact <- fig.contact %>% hide_colorbar() 
fig.contact <- fig.contact %>% colorbar(title = "Vaccine \nEffectiveness (1-RR)")
fig.contact <- fig.contact%>% layout(xaxis = list(title = "VE against susceptibility"), 
                                     yaxis = list(title = "% Contact Increases of \nVaccinated with Vaccinated")) 

fig.contact

#Save as .eps
#orca(fig.contact, "contact_contour", format="eps", width=1200)


######SUPPLEMENTARY MATERIAL##########
#Code for Web Figure 2,3, and 5 and Web Appendix 3
#Code for Web Figure 4 can be found in contact_and_NegVE_SIR.R and contour_data_generation_SIR.R

###DATA GENERATION FOR WEB FIGURE 2 

#VE function for baseline and higher contact with longer timespan
#THREE INPUTS:
#vaccine.type is either "AllOrNothing" or "Leaky", default is "AllOrNothing"
#contact.type is either "Homog." or "Heterog.", default is "Homog." 
#contact.increase is how much higher contact is between vaccinated with vaccinated (i.e. level of vaccinated contact heterogeneity) 
#OUTPUT: Dataset of VE measurements and prop infected from an SEIR model using 300 timesteps
#Dataset contains simulations for VE_S (vaccine efficacy against susc.) and VE_I (vaccine efficacy against infect.) at 0.7 and 0.9


VE_function_highVE<-function(vaccine.type="AllOrNothing", contact.type="Homog.", contact.increase=1.5){
  
  
  names<-c("time","S1","S2","E1","E2","I1","I2","R1","R2",
           "CumIncid1","CumIncid2","Incid1","Incid2",
           "a","b","c", "d","OR.overall",
           "RR.overall","prop.unvac.overall","prop.vac.overall", "prop.total.overall",
           "e","f","g","h", "OR.perday","RR.perday","prop.unvac.perday","prop.vac.perday",
           "label.ve.susc","label.ve.infect")
  
  #Empty df for all VE combos
  Main.VE.df<-data.frame(matrix(ncol=length(names),nrow=0))
  
  #Add column names to df
  colnames(Main.VE.df)<-names
  
  #Full ODE
  ve_model<-function(t,x,parms){ 
    
    
    ncompartment <- 6 # 6 compartments (S, E, I, R)
    ngroup <- length(x)/ncompartment   #number of groups
    S      <- as.matrix(x[1:ngroup])                          #compartment 1
    E    <- as.matrix(x[(ngroup+1):(2*ngroup)])               #compartment 2
    I    <- as.matrix(x[(2*ngroup+1):(3*ngroup)])             #compartment 3
    R   <- as.matrix(x[(3*ngroup+1):(4*ngroup)])              #compartment 4
    CumIncid    <- as.matrix(x[(4*ngroup+1):(5*ngroup)])      #compartment 5
    Incid   <- as.matrix(x[(5*ngroup+1):(6*ngroup)])          #compartment 6
    
    with(as.list(parms),{
      
      N <- S+E+I+R #this just creates a proportional metric per group    
      rho=c(beta,beta*alpha) 
      sigma.matrix=matrix(c(1,1,sigma,sigma),nrow=2,ncol=2)
      
      
      # ODEs
      
      dS <- -as.matrix(S)*rho*(as.matrix(C)*sigma.matrix)%*%as.matrix(I/N)
      
      dE <- as.matrix(S)*rho*(as.matrix(C)*sigma.matrix)%*%as.matrix(I/N) - mu*as.matrix(E) #we assume exposed state isn't infectious
      
      dI <- +mu*as.matrix(E) - gamma*as.matrix(I)       
      
      dR <- +gamma*as.matrix(I)
      
      dCumIncid <-+as.matrix(S)*rho*(as.matrix(C)*sigma.matrix)%*%as.matrix(I/N)
      
      dIncid <-+as.matrix(S)*rho*(as.matrix(C)*sigma.matrix)%*%as.matrix(I/N)  - as.matrix(Incid)
      
      
      #the output has to be in the same order as the model compartments in the begining of the function
      out = c(dS,
              dE,
              dI,         
              dR,          
              dCumIncid,
              dIncid)
      list(out)
      
    }) 
  }
  
  ###General Model Set-up (Applicable to all Model Runs)
  #SET-UP MODEL
  # initial conditions
  npop   <- 100000            # 100,000 total population 
  f      <- c(0.25, 0.75)     # 2 subgroups f[1]=  unvaccinated; f[2] = vaccinated
  Ntot   <- npop*f            # population size in each subgroup (vector of 2 elements)
  ngroup <- length(f)         # length of vector representing number of subgroups
  I_init    <- c(1,1)         # 1 infectious individual per group
  E_init    <- c(0,0)         # 0 exposed individuals
  
  
  cr  = c(6,6)       #baseline number of contacts per day among unvac[1] and vack[2]
  epsilon = 1      #where 1 = proportionate; and 0 = assortative (who contacts whom)
  
  rho_unvac_unvac = (1-epsilon) + epsilon*(Ntot[1]*cr[1]/(Ntot[1]*cr[1] + Ntot[2]*cr[2]))
  rho_unvac_vac  =               epsilon*(Ntot[2]*cr[2]/(Ntot[1]*cr[1] + Ntot[2]*cr[2]))
  rho_vac_unvac   =             epsilon*(Ntot[1]*cr[1]/(Ntot[1]*cr[1] + Ntot[2]*cr[2]))
  rho_vac_vac  =(1-epsilon) + epsilon*(Ntot[2]*cr[2]/(Ntot[1]*cr[1] + Ntot[2]*cr[2]))
  
  weighted_average_cr = (cr[1]*Ntot[1] + cr[2]*Ntot[2])/(Ntot[1]+Ntot[2])
  
  
  #Homogenius or Heterogeneous
  C   <- matrix(0,nrow=ngroup,ncol=ngroup) # contact/mixing matrix of ngroup by ngroup
  
  contact.increase<-contact.increase #set in function
  
  if (contact.type=="Homog."){
    C[1,1]  <- cr[1]*rho_unvac_unvac  # contacts between unvac and unvac
    C[1,2]  <- cr[1]*rho_unvac_vac   # contacts between unvac and vac
    C[2,1]  <- cr[2]*rho_vac_unvac   # contacts between vac and unvac
    C[2,2]  <- cr[2]*rho_vac_vac    # contacts between vac and vac
  }else if (contact.type=="Heterog."){
    C[1,1]  <- cr[1]*rho_unvac_unvac # contacts between unvac and unvac
    C[1,2]  <- cr[1]*rho_unvac_vac   # contacts between unvac and vac
    C[2,1]  <- cr[2]*rho_vac_unvac   # contacts between vac and unvac
    C[2,2]  <- cr[2]*rho_vac_vac*contact.increase   # contacts between vac and vac
  }
  
  #specify other parameters
  gamma   <- 1/10 #recovery rate
  mu<- 1/4 #latency period #2 to 4 days (or 3-4)
  beta<-0.1 #this is for an R0 of 6 
  
  
  #VACCINE EFFICACY ITERATIONS
  ve.start<-0.7
  ve.susc<-ve.start
  
  #Vaccine Efficacy Against Infectiousness Loop
  for (i in 1:2){
    
    
    ve.infect<-ve.start
    
    #Vaccine Efficacy Against Susceptibility Loop
    for(j in 1:2){
      
      
      sigma<-1-ve.infect 
      
      #Here we incorporate VE against suscept. depending on Vaccine Type (all or nothing or leaky)
      if(vaccine.type=="AllOrNothing"){
        alpha=(1-0) 
        R_init<- c(0,Ntot[2]*ve.susc)  
        
      }else if(vaccine.type=="Leaky"){
        alpha=(1-ve.susc)
        R_init<- c(0,0)  
      }
      
      #ODE INitial Conditions & Pars
      S_init    <- Ntot - E_init - I_init - R_init   # initial number susceptible, subtracting number of exposed, infectious, and immune from in each group
      
      rho=c(beta,beta*alpha) 
      sigma.matrix=matrix(c(1,1,sigma,sigma),nrow=2,ncol=2)
      
      CumIncid_init <- I_init
      Incid_init <- I_init
      
      parm <-c(gamma=gamma,
               mu=mu,
               beta=beta, 
               alpha=alpha,
               sigma=sigma,
               C=C)
      
      init <-c(S = S_init,
               E = E_init,
               I = I_init,           
               R = R_init,               
               CumIncid = CumIncid_init,
               Incid = Incid_init)  
      
      
      time = seq(0,300,1)  #run model for 175 days
      
      #RUN MODEL
      
      model.output<-as.data.frame(lsoda(y=init,times=time,func=ve_model,parms=parm)) 
      
      
      
      #Calculating Cumulative Incidence OR and RR
      model.output$a<-model.output$CumIncid2
      model.output$b<-Ntot[2] -model.output$CumIncid2 
      model.output$c<-model.output$CumIncid1
      model.output$d<-Ntot[1] -model.output$CumIncid1
      
      model.output$OR.overall<-1-((model.output$a*model.output$d)/(model.output$b*model.output$c))
      #model.output$RR.overall<-1-((model.output$a/(model.output$a + model.output$b))/(model.output$c/(model.output$c + model.output$d)))
      model.output$RR.overall<-1-((model.output$a/Ntot[2])/(model.output$c/Ntot[1]))
      
      #Calculating Cumulative Incidence Proportion Infected
      model.output$prop.unvac.overall<-model.output$CumIncid1/Ntot[1]
      model.output$prop.vac.overall<-model.output$CumIncid2/Ntot[2]
      model.output$prop.total.overall<-(model.output$CumIncid1 + model.output$CumIncid2)/sum(Ntot)
      
      #Calculating Incidence OR and RR
      model.output$e<-model.output$Incid2
      model.output$f<-Ntot[2] -model.output$Incid2
      model.output$g<-model.output$Incid1
      model.output$h<-Ntot[1] -model.output$Incid1
      
      model.output$OR.perday<-1 - ((model.output$e*model.output$h)/(model.output$f*model.output$g))
      model.output$RR.perday<-1-((model.output$e/(model.output$e + model.output$f))/(model.output$g/(model.output$g + model.output$h)))
      
      #Calculating Incidence Proportion Infected
      model.output$prop.unvac.perday<-model.output$I1/Ntot[1]
      model.output$prop.vac.perday<-model.output$I2/Ntot[2]
      
      #create variable indicator for VEs
      model.output$label.ve.susc=paste("VE Susc.=", ve.susc,sep=" ")
      model.output$label.ve.infect=paste("VE Infect.=", ve.infect, sep=" ")
      
      #Add results to main dataset
      Main.VE.df<-merge(Main.VE.df,model.output,all.x=T,all.y=T,sort=F)
      
      ve.infect<-ve.infect + 0.2
    }
    ve.susc<-ve.susc + 0.2
    
  }
  
  #Add to Vaccine Type and whether contact is Heterog or Homog.
  Main.VE.df$vaccine.type<-paste(vaccine.type)
  Main.VE.df$contact.type<-paste(contact.type)
  
  
  return(Main.VE.df)
  
}


##Create two contact scenario datasets using VE_function_highVE (using 0.7 and 0.9 levels of VES and VEI)

#Baseline contact dataset (50% higher vaccinated contact heterogeneity)
AoN_Homog_VE_highVE_regcontact<-VE_function_highVE(vaccine.type="AllOrNothing", contact.type="Homog.", contact.increase=1.5)
AoN_Heterog_VE_highVE_regcontact<-VE_function_highVE(vaccine.type="AllOrNothing", contact.type="Heterog.", contact.increase=1.5)

subset.homhet_regcontact<-rbind(AoN_Homog_VE_highVE_regcontact,AoN_Heterog_VE_highVE_regcontact)

#Higher Contact dataset (100% higher vaccinated contact heterogeneity)
AoN_Homog_VE_highVE_highcontact<-VE_function_highVE(vaccine.type="AllOrNothing", contact.type="Homog.", contact.increase=2)
AoN_Heterog_VE_highVE_highcontact<-VE_function_highVE(vaccine.type="AllOrNothing", contact.type="Heterog.", contact.increase=2)

subset.homhet_highcontact<-rbind(AoN_Homog_VE_highVE_highcontact,AoN_Heterog_VE_highVE_highcontact)

#Create label groupings
subset.homhet_regcontact$indicator.contact.susc<-factor(paste(subset.homhet_regcontact$contact.type, subset.homhet_regcontact$label.ve.susc, sep="; "), 
                                                        levels=c("Homog.; VE Susc.= 0.7","Homog.; VE Susc.= 0.9", "Heterog.; VE Susc.= 0.7", "Heterog.; VE Susc.= 0.9"))
subset.homhet_highcontact$indicator.contact.susc<-factor(paste(subset.homhet_highcontact$contact.type, subset.homhet_highcontact$label.ve.susc, sep="; "),
                                                         levels=c("Homog.; VE Susc.= 0.7","Homog.; VE Susc.= 0.9", "Heterog.; VE Susc.= 0.7", "Heterog.; VE Susc.= 0.9"))

##Create Proportion Infected variable for both high and low contact scenarios
subset.homhet_regcontact$prop.infected<-(subset.homhet_regcontact$E1 + subset.homhet_regcontact$E2 + subset.homhet_regcontact$I1 + subset.homhet_regcontact$I2)/
  (subset.homhet_regcontact$E1 + subset.homhet_regcontact$E2 + subset.homhet_regcontact$I1 + subset.homhet_regcontact$I2 + subset.homhet_regcontact$S1 + subset.homhet_regcontact$S2 + subset.homhet_regcontact$R1 + subset.homhet_regcontact$R2)

subset.homhet_highcontact$prop.infected<-(subset.homhet_highcontact$E1 + subset.homhet_highcontact$E2 + subset.homhet_highcontact$I1 + subset.homhet_highcontact$I2)/
  (subset.homhet_highcontact$E1 + subset.homhet_highcontact$E2 + subset.homhet_highcontact$I1 + subset.homhet_highcontact$I2 + subset.homhet_highcontact$S1 + subset.homhet_highcontact$S2 + subset.homhet_highcontact$R1 + subset.homhet_highcontact$R2)


### SUPPLEMENTARY WEB FIGURE 2 LINE PLOTS (VE using Relative Risk and Proportion Infected for High VE Scenarios)

colours.hom.het<-c( "orangered1", "darkred","slateblue1", "navyblue")
lines.hom.het<-c("solid", "twodash")

##VE and Prop Infected using baseline vaccinated contact heterogeneity (50% higher)

#Vaccine Effectiveness (1-RR) Plot
Ve.RR.overall.plot_regcontact<-ggplot(subset.homhet_regcontact, aes(x=time, y=RR.overall, color=indicator.contact.susc, linetype=label.ve.infect)) +
  theme_classic()  +  geom_hline(yintercept=0, linetype="dotted", color = "darkgrey", size=0.8) +
  geom_line(size=1.2) + ggtitle("Vaccine Effectiveness and Proportion Infected Across Time") + 
  labs(y="Vaccine Effectiveness (1-Relative Risk)", x = "Days", color="VE Values") + 
  scale_x_continuous(limits = c(0, 300),breaks=seq(0,300,50) ) + scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values=colours.hom.het) + scale_linetype_manual(values=lines.hom.het) +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#Prop. Infected
prop.infected.plot_regcontact<-ggplot(subset.homhet_regcontact, aes(x=time, y=prop.infected, color=indicator.contact.susc, linetype=label.ve.infect)) +
  theme_classic() +
  geom_line(size=1.2) +
  labs(y="Proportion Infected", x = "Days", color="VE Values") + 
  scale_x_continuous(limits = c(0, 300), breaks=seq(0,300,50) ) + 
  scale_y_continuous(limits = c(0, 0.5), breaks=seq(0,0.5,0.1) )+ 
  scale_color_manual(values=colours.hom.het) + scale_linetype_manual(values=lines.hom.het)

##VE and Prop Infected using high vaccinated contact heterogeneity (100% higher vaccinated contact heterogeneity)

#Vaccine Effectiveness (1-RR) Plot 
Ve.RR.overall.plot_highcontact<-ggplot(subset.homhet_highcontact, aes(x=time, y=RR.overall, color=indicator.contact.susc, linetype=label.ve.infect)) +
  theme_classic()  +  geom_hline(yintercept=0, linetype="dotted", color = "darkgrey", size=0.8) +
  geom_line(size=1.2) + ggtitle("Vaccine Effectiveness and Proportion Infected Across Time") + 
  labs(y="Vaccine Effectiveness (1-Relative Risk)", x = "Days", color="VE Values") + 
  scale_x_continuous(limits = c(0, 300),breaks=seq(0,300,50) ) + scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values=colours.hom.het) + scale_linetype_manual(values=lines.hom.het) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),
        axis.ticks.x=element_blank(), axis.ticks.y=element_blank())

#Prop Infected
prop.infected.plot_highcontact<-ggplot(subset.homhet_highcontact, aes(x=time, y=prop.infected, color=indicator.contact.susc, linetype=label.ve.infect)) +
  theme_classic() +
  geom_line(size=1.2) +
  labs(y="Proportion Infected", x = "Days", color="VE Values") + 
  scale_x_continuous(limits = c(0, 300), breaks=seq(0,300,50) ) + 
  scale_y_continuous(limits = c(0, 0.5), breaks=seq(0,0.5,0.1) )+ 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + scale_color_manual(values=colours.hom.het) + scale_linetype_manual(values=lines.hom.het)




line.plots.overall.highVE<-ggarrange(Ve.RR.overall.plot_regcontact, Ve.RR.overall.plot_highcontact,
                                     prop.infected.plot_regcontact,prop.infected.plot_highcontact,
                                     ncol=2, nrow=2, common.legend=T, legend = "right")

line.plots.overall.highVE


### INFORMATION FOR SUPPLEMENTARY WEB APPENDIX 3
#Exploring Vaccine Effectiveness crossover for Figure 1 (Heterogeneous Contact with low VEs)

#Data from VE_function above
ve.all01<-subset(subset.homhet,indicator.contact.susc=="Heterog.; VE Susc.= 0.1" & label.ve.infect == "VE Infect.= 0.1" )
ve.01.05<-subset(subset.homhet,indicator.contact.susc=="Heterog.; VE Susc.= 0.1" & label.ve.infect == "VE Infect.= 0.5")

#Reminder of our pop size and vaccination prop
npop   <- 100000            # 100,000 total population 
f      <- c(0.25, 0.75)     # 2 subgroups f[1]=  unvaccinated; f[2] = vaccinated
Ntot   <- npop*f 

#calculate proportion of susc. unvac and vac for two scenarios
#VES=0.1 and VEI=0.1
ve.all01$prop.susc.S1<-ve.all01$S1/Ntot[1]
ve.all01$prop.susc.S2<-ve.all01$S2/Ntot[2]

#VES=0.1 and VEI=0.5
ve.01.05$prop.susc.S1<-ve.01.05$S1/Ntot[1]
ve.01.05$prop.susc.S2<-ve.01.05$S2/Ntot[2]

#Add VES=0.1 to proportion of susceptible vaccinated
ve.all01$prop.susc.S2.alpha<-ve.all01$prop.susc.S2 + 0.1
ve.01.05$prop.susc.S2.alpha<-ve.01.05$prop.susc.S2 + 0.1

#calculate the difference between prop. susc. vac. + VES and prop. susc. unvac.
ve.all01$susc.diff<-  ve.all01$prop.susc.S2.alpha - ve.all01$prop.susc.S1
ve.01.05$susc.diff<-  ve.01.05$prop.susc.S2.alpha - ve.01.05$prop.susc.S1

#For Figure 1 - calculating the crossovers
ve.all01[c("time","susc.diff", "RR.overall")] #VE becomes positive at t=51 (with VES=VEI=0.1)
ve.01.05[c("time","susc.diff", "RR.overall")] #VE becomes positive at t=75 (with VES=0.1 and VEI=0.5)

ve.all01[c("time","prop.susc.S2", "prop.susc.S1", "RR.overall")]

#The underestimate is maximized in the epidemic growth state
min(ve.all01["RR.overall"])
ve.all01$time[ve.all01$RR.overall==min(ve.all01["RR.overall"])]
min(ve.01.05["RR.overall"])
ve.01.05$time[ve.01.05$RR.overall==min(ve.01.05["RR.overall"])]


###SUPPLEMENTARY WEB FIGURE 3 DIFFERENCE IN PROP SUSCPT VAC WITH VES and PROP SUSCP UNVAC

#Main Figure
susc.diff.plot<- ggplot(ve.01.05, aes(x=time, y=susc.diff)) +
  theme_classic() +
  geom_hline(yintercept=0, linetype="dotted", color = "darkgrey", size=0.8) +
  geom_line(size=0.9, color="#56B4E9", linetype="twodash") +
  geom_line(data=ve.all01,aes(y = susc.diff), color = "#56B4E9", linetype = "solid",size=0.9) +
  scale_x_continuous(limits=c(0,175),breaks=seq(0,175,25) ) +  scale_y_continuous(limits = c(-0.06, 0.10),breaks=seq(-0.08,0.12,0.02))+ 
  labs(y="Difference in Prop. Sucept. with VE[S] Between Vacc. and Unvacc. Groups", x = "Time") 

susc.diff.plot

#Zoomed in plot:
zoomed.susc.diff.plot<-ggplot(ve.01.05, aes(x=time, y=susc.diff)) +
  theme_classic() +
  geom_hline(yintercept=0, linetype="dotted", color = "darkgrey", size=0.8) +
  geom_line(size=0.9, color="#56B4E9", linetype="twodash") +
  geom_line(data=ve.all01,aes(y = susc.diff), color = "#56B4E9", linetype = "solid",size=0.9) +
  scale_x_continuous(limits=c(0,15),breaks=seq(0,15,5) ) +  scale_y_continuous(limits = c(-0.0002, 0.0002),breaks=seq(-0.0002,0.0002,0.0002))+ 
  labs(y="Difference in Prop. Sucept. with VE[S] Between Vacc. and Unvacc. Groups", x = "Time") 

zoomed.susc.diff.plot #Warning appears as only taking a subset of the data

###SUPPLEMENTARY WEB FIGURE 4 RESULTS FOUND IN SIR CODING SCRIPTS 
#SIR CODE FILE NAMES: contact_and_NegVE_SIR.R and contour_data_generation_SIR.R

###SUPPLEMENTARY WEB FIGURE 5 - UNDERESTIMATION CONTOUR PLOTS

#Read in contour plot data from contour_data_generation_SEIR.R
contour.under_0.1<-read.csv("data/VEI_and_ContactIncrease_Under_VES0_1.csv")
contour.under_0.3<-read.csv("data/VEI_and_ContactIncrease_Under_VES0_3.csv")
contour.under_0.5<-read.csv("data/VEI_and_ContactIncrease_Under_VES0_5.csv")
contour.under_0.7<-read.csv("data/VEI_and_ContactIncrease_Under_VES0_7.csv")
contour.under_0.9<-read.csv("data/VEI_and_ContactIncrease_Under_VES0_9.csv")

#calculate the maximum degree of underestimation
contour.under_0.1$max.underestimate<-0.1 - contour.under_0.1$min.RR.overall
contour.under_0.3$max.underestimate<-0.3 - contour.under_0.3$min.RR.overall
contour.under_0.5$max.underestimate<-0.5 - contour.under_0.5$min.RR.overall
contour.under_0.7$max.underestimate<-0.7 - contour.under_0.7$min.RR.overall
contour.under_0.9$max.underestimate<-0.9 - contour.under_0.9$min.RR.overall

#create contact increase % 
contour.under_0.1$contact.increase.percent<-(contour.under_0.1$contact.increase - 1)*100
contour.under_0.3$contact.increase.percent<-(contour.under_0.3$contact.increase - 1)*100
contour.under_0.5$contact.increase.percent<-(contour.under_0.5$contact.increase - 1)*100
contour.under_0.7$contact.increase.percent<-(contour.under_0.7$contact.increase - 1)*100
contour.under_0.9$contact.increase.percent<-(contour.under_0.9$contact.increase - 1)*100

#Maximum Underestimation contour plot (VES=0.1)
fig.contour.under_0.1<-plot_ly(z = contour.under_0.1$max.underestimate,x=contour.under_0.1$ve.infect, 
                               y=contour.under_0.1$contact.increase.percent,
                               type = "contour", autocontour=F, line = list(width = 1, color = "black"),
                               contours = list(coloring = 'heatmap',  start = 0, end = 0.75, size = 0.05, showlabels=T,
                                               labelfont = list(size = 27)),
                               colorscale='Hot', reversescale=T,
                               line = list(smoothing = 0.3))
fig.contour.under_0.1 <- fig.contour.under_0.1 %>% colorbar(title = "Max. Underestimate \nof Vaccine \nEffectiveness")
fig.contour.under_0.1 <- fig.contour.under_0.1 %>% layout(xaxis = list(title = "VE against infectiousness"), 
                                                          yaxis = list(title = "% Contact Increases of \nVaccinated with Vaccinated")) 

fig.contour.under_0.1

#Save as .eps
#orca(fig.contour.under_0.1, "fig.under_01", format="eps", width=800)

#Maximum Underestimation contour plot (VES=0.3)
fig.contour.under_0.3 <-plot_ly(z = contour.under_0.3$max.underestimate,x=contour.under_0.3$ve.infect, 
                                y=contour.under_0.3$contact.increase.percent,
                                type = "contour", autocontour=F, line = list(width = 1, color = "black"),
                                contours = list(coloring = 'heatmap',  start = 0, end = 0.75, size = 0.05, showlabels=T,
                                                labelfont = list(size = 27)),
                                colorscale="Hot", reversescale=T,
                                line = list(smoothing = 0.3))
fig.contour.under_0.3 <- fig.contour.under_0.3 %>% colorbar(title = "Max. Underestimate \nof Vaccine \nEffectiveness")
fig.contour.under_0.3 <- fig.contour.under_0.3 %>% layout(xaxis = list(title = "VE against infectiousness"), 
                                                          yaxis = list(title = "% Contact Increases of \nVaccinated with Vaccinated")) 


fig.contour.under_0.3
#Save as .eps
#orca(fig.contour.under_0.3, "fig.under_03", format="eps", width=800)

#Maximum Underestimation contour plot (VES=0.5)
fig.contour.under_0.5 <-plot_ly(z = contour.under_0.5$max.underestimate,x=contour.under_0.5$ve.infect, 
                                y=contour.under_0.5$contact.increase.percent,
                                type = "contour", autocontour=F, line = list(width = 1, color = "black"),
                                contours = list(coloring = 'heatmap',  start = 0, end = 0.75, size = 0.05, showlabels=T,
                                                labelfont = list(size = 27)),
                                colorscale="Hot", reversescale=T,
                                line = list(smoothing = 0.3))
fig.contour.under_0.5 <- fig.contour.under_0.5 %>% colorbar(title =  "Max. Underestimate \nof Vaccine \nEffectiveness")
fig.contour.under_0.5 <- fig.contour.under_0.5 %>% layout(xaxis = list(title = "VE against infectiousness"), 
                                                          yaxis = list(title = "% Contact Increases of \nVaccinated with Vaccinated")) 


fig.contour.under_0.5

#Save as .eps
#orca(fig.contour.under_0.5, "fig.under_05", format="eps", width=800)


#Maximum Underestimation contour plot (VES=0.7)
fig.contour.under_0.7 <-plot_ly(z = contour.under_0.7$max.underestimate,x=contour.under_0.7$ve.infect, 
                                y=contour.under_0.7$contact.increase.percent,
                                type = "contour", autocontour=F, line = list(width = 1, color = "black"),
                                contours = list(coloring = 'heatmap',  start = 0, end = 0.75, size = 0.05, showlabels=T,
                                                labelfont = list(size = 27)),
                                colorscale="Hot", reversescale=T,
                                line = list(smoothing = 0.3))
fig.contour.under_0.7 <- fig.contour.under_0.7 %>% colorbar(title =  "Max. Underestimate \nof Vaccine \nEffectiveness")
fig.contour.under_0.7 <- fig.contour.under_0.7 %>% layout(xaxis = list(title = "VE against infectiousness"), 
                                                          yaxis = list(title = "% Contact Increases of \nVaccinated with Vaccinated")) 


fig.contour.under_0.7
#Save as .eps 
#orca(fig.contour.under_0.7, "fig.under_07", format="eps", width=800)


#Maximum Underestimation contour plot (VES=0.9)
fig.contour.under_0.9 <-plot_ly(z = contour.under_0.9$max.underestimate,x=contour.under_0.9$ve.infect, 
                                y=contour.under_0.9$contact.increase.percent,
                                type = "contour", autocontour=F, line = list(width = 1, color = "black"),
                                contours = list(coloring = 'heatmap',  start = 0, end = 0.75, size = 0.05, showlabels=T,
                                                labelfont = list(size = 27)),
                                colorscale="Hot", reversescale=T,
                                line = list(smoothing = 0.3))
fig.contour.under_0.9 <- fig.contour.under_0.9 %>% colorbar(title =  "Max. Underestimate \nof Vaccine \nEffectiveness")
fig.contour.under_0.9 <- fig.contour.under_0.9 %>% layout(xaxis = list(title = "VE against infectiousness"), 
                                                          yaxis = list(title = "% Contact Increases of \nVaccinated with Vaccinated")) 


fig.contour.under_0.9

#Save as .eps
#orca(fig.contour.under_0.9, "fig.under_09", format="eps", width=800)

