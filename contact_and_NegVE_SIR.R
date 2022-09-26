###################################################################### #
#                                                                      #                                                   #
# Vaccinated Contact Heterogeneity and Negative Vaccine Effectiveness  #
#                 SUPP MAT CODE: SIR Model Results                     #
#                                                                      #
########################################################################

#PACKAGES
require(deSolve)
require(tidyverse)
require(ggplot2)
require(reshape)
library(ggpubr)
library(plotly)

#CODE FOR SUPPLEMENTARY FIGURE 4 (SIR MODEL RESULTS)

#VE FUNCTION HAS TWO INPUTS:
#vaccine.type is either "AllOrNothing" or "Leaky", default is "AllOrNothing"
#contact.type is either "Homog." (random mixing) or "Heterog." (vaccinated contact heterogeneity)
  #default is "Homog.
#OUTPUTS: dataset with VE measurements and proportion infected using an SIR model
#EXTRA DETAILS: 
#vaccine efficacy against susceptibility (VE_S) and vaccine efficacy against infectiousness (VE_I) are 0.1, 0.5, 0.9
#vaccinated contact heterogeneity assumed to be 50% higher than random mixing case

VE_function<-function(vaccine.type="AllOrNothing", contact.type="Homog."){
  
  #Empty df for all VE combos
  Main.VE.df<-data.frame(matrix(ncol=40,nrow=0))
  
  names<-c("time","S1","S2","I1","I2","R1","R2","CumIncid1","CumIncid2","Incid1","Incid2",
           "CumFOI1","CumFOI2","FOI1","FOI2", "P.time.infected1","P.time.infected2","Reff1",
           "Reff2","Reff.avg","FOI.avg","a","b","c", "d","OR.overall",
           "RR.overall","prop.unvac.overall","prop.vac.overall", "prop.total.overall",
           "e","f","g","h", "OR.perday","RR.perday","prop.unvac.perday","prop.vac.perday",
           "label.ve.susc","label.ve.infect")
  colnames(Main.VE.df)<-names
  
  #Full ODE
  ve_model<-function(t,x,parms){ 
    
    
    ncompartment <- 9 # 5 compartments (S, E, I, R)
    ngroup <- length(x)/ncompartment   #number of groups
    S      <- as.matrix(x[1:ngroup])                          #compartment 1
    I    <- as.matrix(x[(ngroup+1):(2*ngroup)])               #compartment 2
    R   <- as.matrix(x[(2*ngroup+1):(3*ngroup)])              #compartment 3
    CumIncid    <- as.matrix(x[(3*ngroup+1):(4*ngroup)])      #compartment 4
    Incid   <- as.matrix(x[(4*ngroup+1):(5*ngroup)])          #compartment 5
    CumFOI  <- as.matrix(x[(5*ngroup+1):(6*ngroup)])          #compartment 6
    FOI  <- as.matrix(x[(6*ngroup+1):(7*ngroup)])             #compartment 7
    P.time.infected  <- as.matrix(x[(7*ngroup+1):(8*ngroup)]) #compartment 8
    Reff  <- as.matrix(x[(8*ngroup+1):(9*ngroup)])            #compartment 9
    
    with(as.list(parms),{
      
      N <- S+I+R #this just creates a proportional metric per group    
      rho=c(beta,beta*alpha) 
      sigma.matrix=matrix(c(1,1,sigma,sigma),nrow=2,ncol=2)
      
      
      # ODEs
      
      dS <- -as.matrix(S)*rho*(as.matrix(C)*sigma.matrix)%*%as.matrix(I/N)
      
      dI <- as.matrix(S)*rho*(as.matrix(C)*sigma.matrix)%*%as.matrix(I/N) - gamma*as.matrix(I)       
      
      dR <- +gamma*as.matrix(I)
      
      dCumIncid <-+as.matrix(S)*rho*(as.matrix(C)*sigma.matrix)%*%as.matrix(I/N)
      
      dIncid <-+as.matrix(S)*rho*(as.matrix(C)*sigma.matrix)%*%as.matrix(I/N)  - as.matrix(Incid)
      
      dCumFOI <-+rho*(as.matrix(C)*sigma.matrix)%*%as.matrix(I/N) 
      
      dFOI <-+rho*(as.matrix(C)*sigma.matrix)%*%as.matrix(I/N) - as.matrix(FOI)
      
      dP.time.infected<-t*as.matrix(S)*rho*(as.matrix(C)*sigma.matrix)%*%as.matrix(I/N)
      
      dReff<-as.matrix(S/N)*rho*(as.matrix(C)*sigma.matrix)%*%as.matrix(c(1,1))*1/gamma - as.matrix(Reff)
      
      #the output has to be in the same order as the model compartments in the begining of the function
      out = c(dS,           
              dI,         
              dR,          
              dCumIncid,
              dIncid,
              dCumFOI,
              dFOI,
              dP.time.infected,
              dReff)
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
  I_init    <- c(1,1)  # 1 wt infectious individual
  
  
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
  beta<-0.1 #this is for an R0 of 6 
  
  
  #VACCINE EFFICACY ITERATIONS
  ve.start<-0.1
  ve.susc<-ve.start
  
  #Vaccine Efficacy Against Infectiousness Loop
  for (i in 1:3){
    
    
    ve.infect<-ve.start
    
    #Vaccine Efficacy Against Susceptibility Loop
    for(j in 1:3){
      
      
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
      S_init    <- Ntot - I_init - R_init   # initial number susceptible, subtracting number infectious in each group
      CumFOI_init <- c(beta,beta*alpha)*C %*% as.matrix(I_init/Ntot) ####double check this
      FOI_init <- c(beta,beta*alpha)*C %*% as.matrix(I_init/Ntot)
      
      rho=c(beta,beta*alpha) 
      sigma.matrix=matrix(c(1,1,sigma,sigma),nrow=2,ncol=2)
      Reff_init<-as.matrix(S_init/Ntot)*rho*(as.matrix(C)*sigma.matrix)%*%as.matrix(c(1,1))*1/gamma
      
      CumIncid_init <- I_init
      Incid_init <- I_init
      P.time.infected_init<-rep(0,ngroup)
      
      parm <-c(gamma=gamma,
               beta=beta, 
               alpha=alpha,
               sigma=sigma,
               C=C)
      
      init <-c(S = S_init,               
               I = I_init,           
               R = R_init,               
               CumIncid = CumIncid_init,
               Incid = Incid_init,
               CumFOI = CumFOI_init,
               FOI = FOI_init,
               P.time.infected=P.time.infected_init,
               Reff=Reff_init)  
      
      
      time = seq(0,150,1)  #run model for 150 days
      
      #RUN MODEL
      
      model.output<-as.data.frame(lsoda(y=init,times=time,func=ve_model,parms=parm)) 
      
      
      
      #Calculating Reffective and FOI average
      model.output$Reff.avg<-model.output$Reff1*f[1] + model.output$Reff2*f[2]
      model.output$FOI.avg<-model.output$FOI1*f[1] + model.output$FOI2*f[2]
      
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


#Create Datasets for VE_function
AoN_Homog_VE<-VE_function(vaccine.type="AllOrNothing", contact.type="Homog.")
AoN_Heterog_VE<-VE_function(vaccine.type="AllOrNothing", contact.type="Heterog.")


#Create a single dataset with both homog. and heterog together (remove VE = 0.9)
subset.Homog<-subset(AoN_Homog_VE, label.ve.infect !="VE Infect.= 0.9" & label.ve.susc !="VE Susc.= 0.9")
subset.Heterog<-subset(AoN_Heterog_VE, label.ve.infect !="VE Infect.= 0.9" & label.ve.susc !="VE Susc.= 0.9")

subset.homhet<-rbind(subset.Homog, subset.Heterog)

#Create label grouping Homog with Susc
subset.homhet$indicator.contact.susc<-paste(subset.homhet$contact.type, subset.homhet$label.ve.susc, sep="; ")

subset.homhet$prop.infected<-(subset.homhet$I1 + subset.homhet$I2)/(subset.homhet$I1 + subset.homhet$I2 + subset.homhet$S1 + subset.homhet$S2 + subset.homhet$R1 + subset.homhet$R2)


###FIGURE 1 (a and b) LINE PLOTS (Relative Risk and Proportion Infected)

colours.hom.het<-c("#56B4E9", "#0072B2", "#E69F00", "#D55E00")


#Create time and Vaccine Effectiveness (1-RR) Plot
#t=24 and t=39 when RR crosses over in heterog case with low VE
Ve.RR.overall.plot<-ggplot(subset.homhet, aes(x=time, y=RR.overall, color=indicator.contact.susc, linetype=label.ve.infect)) +
  theme_classic() + geom_hline(yintercept=0, linetype="dotted", color = "darkgrey", size=0.8) +
  geom_segment(aes(x=24, y=-Inf, xend=24, yend=0),linetype="dotted", color = "darkgrey", size=0.8) +
  geom_segment(aes(x=39, y=-Inf, xend=39, yend=0),linetype="dotted", color = "darkgrey", size=0.8) +
  geom_line(size=1.2) +
  ggtitle("Vaccine Effectiveness and Proportion Infected Across Time") + 
  labs(y="Vaccine Effectiveness (1-Relative Risk)", x = "Days", color="VE Values") + 
  scale_x_continuous(limits = c(0, 150)) + scale_y_continuous(limits = c(-0.3, 0.8)) +
  scale_color_manual(values=colours.hom.het) + theme(axis.title.x=element_blank(),
                                                     axis.text.x=element_blank(),
                                                     axis.ticks.x=element_blank())

#Create prop.infected
prop.infected.plot<-ggplot(subset.homhet, aes(x=time, y=prop.infected, color=indicator.contact.susc, linetype=label.ve.infect)) +
  theme_classic() +
  geom_line(size=1.2) +
  labs(y="Proportion Infected", x = "Days", color="VE Values") + 
  scale_x_continuous(limits = c(0, 150)) +
  scale_color_manual(values=colours.hom.het) 


line.plots.overall<-ggarrange(Ve.RR.overall.plot, prop.infected.plot, ncol=1, common.legend=T, legend = "right")

line.plots.overall

###FIGURE 1 (c and d) CONTOUR PLOTS (levels of VES with VEI and vaccinated contact heterogeneity)

VeT.contour.df<-read.csv("data/VES_and_VEI_SIR.csv", )
contact.VES.df<-read.csv("data/VES_and_ContactIncrease_SIR.csv")

#Subsetting VES<0.51 as higher values produces the same results as 0.5
veT.susc0.5<-subset(VeT.contour.df, ve.susc<0.51)

#Contour Plot for VES and VEI
fig.veT<-plot_ly(z = veT.susc0.5$min.RR.overall,x=veT.susc0.5$ve.susc, 
                 y=veT.susc0.5$ve.infect, type = "contour", autocontour=F,line = list(width = 1, color = "black"),
                 contours = list(coloring = 'heatmap',  start = -.65, end = 0, size = 0.05, showlabels=T),
                 colorscale='Jet',reversescale=T,
                 line = list(smoothing = 0.3))
fig.veT <- fig.veT %>% colorbar(title = "Vaccine \nEffectiveness (1-RR)")  
fig.veT <- fig.veT%>% layout(xaxis = list(title = "VE against susceptibility"), 
                             yaxis = list(title = "VE against infectiousness")) 


fig.veT


#Subsetting VES<0.51 as higher values produces the same results as 0.5
contact.ve.susc0.5<-subset(contact.VES.df, ve.susc<0.51)

#Convert contact rate into % increase above Homogeneous contact scenario
contact.ve.susc0.5$contact.increase.percent<-(contact.ve.susc0.5$contact.increase - 1)*100

#Contour Plot for VES and vaccinated contact heterogeneity
fig.contact <-plot_ly(z = contact.ve.susc0.5$min.RR.overall,x=contact.ve.susc0.5$ve.susc, 
                      y=contact.ve.susc0.5$contact.increase.percent,
                      type = "contour", autocontour=F, line = list(width = 1, color = "black"),
                      contours = list(coloring = 'heatmap',  start = -.65, end = 0, size = 0.05, showlabels=T),
                      colorscale='Jet',reversescale=T,
                      line = list(smoothing = 0.3))
#fig.contact <- fig.contact %>% hide_colorbar() 
fig.contact <- fig.contact %>% colorbar(title = "Vaccine \nEffectiveness (1-RR)")
fig.contact <- fig.contact%>% layout(xaxis = list(title = "VE against susceptibility"), 
                                     yaxis = list(title = "% Contact Increases of \nVaccinated with Vaccinated")) 

fig.contact


####Exploring Vaccine Effectiveness crossover from negative to positive with vaccinated contact heterogeneity and low VES

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

ve.all01[c("time","susc.diff", "RR.overall")] #VES=VEI=0.1 has a crossover at t=24
ve.01.05[c("time","susc.diff", "RR.overall")] #VES=0.1 and VEI=0.5 has a crossover at t=40

ve.all01[c("time","prop.susc.S2", "prop.susc.S1", "RR.overall")]



#Figure showing crossovers with SIR model (Same general results as SEIR model) 

susc.diff.plot<- ggplot(ve.01.05, aes(x=time, y=susc.diff)) +
  theme_classic() +
  geom_hline(yintercept=0, linetype="dotted", color = "darkgrey", size=0.8) +
  geom_line(size=0.9, color="#56B4E9", linetype="twodash") +
  geom_line(data=ve.all01,aes(y = susc.diff), color = "#56B4E9", linetype = "solid",size=0.9) +
  scale_x_continuous(limits=c(0,150),breaks=seq(0,150,50) ) +  scale_y_continuous(limits = c(-0.06, 0.10),breaks=seq(-0.08,0.12,0.02))+ 
  labs(y="Difference in Prop. Sucept. with VE[S] Between Vacc. and Unvacc. Groups", x = "Time") 

susc.diff.plot

#Zoomed in plot:
zoomed.susc.diff.plot<-ggplot(ve.01.05, aes(x=time, y=susc.diff)) +
  theme_classic() +
  geom_hline(yintercept=0, linetype="dotted", color = "darkgrey", size=0.8) +
  geom_line(size=0.9, color="#56B4E9", linetype="twodash") +
  geom_line(data=ve.all01,aes(y = susc.diff), color = "#56B4E9", linetype = "solid",size=0.9) +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,5) ) +  scale_y_continuous(limits = c(-0.0001, 0.0001),breaks=seq(-0.0001,0.0001,0.0001))+ 
  labs(y="Difference in Prop. Sucept. with VE[S] Between Vacc. and Unvacc. Groups", x = "Time") 

zoomed.susc.diff.plot


