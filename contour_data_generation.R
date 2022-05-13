#Contour function that runs scenarios for different values of VES and VEI
veT.contour<-function(contact.type="Homog.", vaccine.type="AllOrNothing"){
  
  n.iter<-51
  
  Contour.VE.df<-data.frame(matrix(ncol=6,nrow=n.iter*n.iter))
  
  names<-c("min.OR.overall","min.RR.overall",
           "min.OR.perday","min.RR.perday", 
           "ve.susc","ve.infect")
  colnames(Contour.VE.df)<-names
  
  #Full ODE
  neg_ve_model<-function(t,x,parms){ 
    
    
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
      
      N <- S+I+R #this just creates a proportional metric per group  (I think it should be n)     
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
      
      #dReff<-as.matrix(S/N)*rho*(as.matrix(C)*sigma.matrix)%*%as.matrix(c(1,1))*1/gamma - as.matrix(Reff)
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
  
  ######START OF VE
  ve.start=0
  
  ve.susc<-ve.start
  
  ###VE for susceptibility
  for (i in 1:n.iter){
    
    ve.infect<-ve.start
    
    for (j in 1:n.iter){
      
      nrow<-n.iter*(i-1) + j #this is the row number in the main dataset
      
      #SPECIFIC MODEL SET-UP (SET VE)
      
      sigma<-1-ve.infect #sigma is the vaccine efficacy for infectiousness for both ‘leaky’ and ‘all-or-nothing’ vaccines.
      
      #Here we set alpha (vaccine efficacy against suscept.) depending on Vaccine Type
      if(vaccine.type=="AllOrNothing"){
        alpha=(1-0) 
        R_init<- c(0,Ntot[2]*ve.susc)  
        
      }else if(vaccine.type=="Leaky"){
        alpha=(1-ve.susc)
        R_init<- c(0,0)  
      }
      
      S_init    <- Ntot - I_init - R_init   # initial number susceptible, subtracting number infectious in each group
      CumFOI_init <- c(beta,beta*alpha)*C %*% as.matrix(I_init/Ntot) ####double check this
      FOI_init <- c(beta,beta*alpha)*C %*% as.matrix(I_init/Ntot)
      
      rho=c(beta,beta*alpha) 
      sigma.matrix=matrix(c(1,1,sigma,sigma),nrow=2,ncol=2)
      Reff_init<-as.matrix(S_init/Ntot)*rho*(as.matrix(C)*sigma.matrix)%*%as.matrix(c(1,1))*1/gamma
      
      
      CumIncid_init <- I_init
      Incid_init <- I_init
      P.time.infected_init<-rep(0,ngroup)
      
      parms <-c(gamma=gamma,
                beta=beta, 
                alpha=alpha,
                sigma=sigma,
                C=C)
      
      inits <-c(S = S_init,               
                I = I_init,           
                R = R_init,               
                CumIncid = CumIncid_init,
                Incid = Incid_init,
                CumFOI = CumFOI_init,
                FOI = FOI_init,
                P.time.infected=P.time.infected_init,
                Reff=Reff_init)  
      
      
      times = seq(0,300,1)  #run model for 300 days
      
      #RUN MODEL
      model.output= as.data.frame(lsoda(inits,times,neg_ve_model,parms))  
      
      #Calculating Cumulative Incidence OR and RR
      model.output$a<-model.output$CumIncid2
      model.output$b<-Ntot[2] -model.output$CumIncid2 
      model.output$c<-model.output$CumIncid1
      model.output$d<-Ntot[1] -model.output$CumIncid1
      
      model.output$OR.overall<-1-((model.output$a*model.output$d)/(model.output$b*model.output$c))
      model.output$RR.overall<-1-((model.output$a/(model.output$a + model.output$b))/(model.output$c/(model.output$c + model.output$d)))
      
      #Calculating Incidence OR and RR
      model.output$e<-model.output$Incid2
      model.output$f<-Ntot[2] -model.output$Incid2
      model.output$g<-model.output$Incid1
      model.output$h<-Ntot[1] -model.output$Incid1
      
      model.output$OR.perday<-1 - ((model.output$e*model.output$h)/(model.output$f*model.output$g))
      model.output$RR.perday<-1-((model.output$e/(model.output$e + model.output$f))/(model.output$g/(model.output$g + model.output$h)))
      
      #GET MINIMUM VALUES
      
      Contour.VE.df$min.OR.overall[nrow]<-min(model.output$OR.overall)
      Contour.VE.df$min.RR.overall[nrow]<-min(model.output$RR.overall)
      
      Contour.VE.df$min.OR.perday[nrow]<-min(model.output$OR.perday)
      Contour.VE.df$min.RR.perday[nrow]<-min(model.output$RR.perday)
      
      Contour.VE.df$ve.susc[nrow]<-ve.susc
      Contour.VE.df$ve.infect[nrow]<-ve.infect
      
      ve.infect<-ve.infect + 0.02 
    }
    
    ve.susc<-ve.susc + 0.02
  }
  
  return(Contour.VE.df)
}

#Generating scenarios (can take a while to run)
VeT.contour.df<-veT.contour(contact.type="Heterog.", vaccine.type="AllOrNothing")

write.csv(VeT.contour.df,"VES_and_VEI.csv")

#Contour function that runs scenarios for different values of VES and Contact Heterogeneity
veS.contact.contour<-function(vaccine.type="AllOrNothing"){
  
  vaccine.type="AllOrNothing"
  n.iter<-51
  
  
  Contour.VE.df<-data.frame(matrix(ncol=7,nrow=n.iter*n.iter))
  
  names<-c("min.OR.overall","min.RR.overall",
           "min.OR.perday","min.RR.perday", 
           "ve.susc","ve.infect", "contact.increase")
  colnames(Contour.VE.df)<-names
  
  #Full ODE
  ve_model_contact<-function(t,x,parms){ 
    
    
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
      
      N <- S+I+R #this just creates a proportional metric per group  (I think it should be n)     
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
      
      #dReff<-as.matrix(S/N)*rho*(as.matrix(C)*sigma.matrix)%*%as.matrix(c(1,1))*1/gamma - as.matrix(Reff)
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
  
  
  #specify other parameters
  gamma   <- 1/10 #recovery rate
  beta<-0.1 #this is for an R0 of 6 
  
  ######VE
  ve.start=0
  
  ve.susc<-ve.start
  ve.infect<-0
  
  ###VE for susceptibility
  for (i in 1:n.iter){
    
    contact.increase<-1
    
    for (j in 1:n.iter){
      
      nrow<-n.iter*(i-1) + j #this is the row number in the main dataset
      
      #Contact Matrix
      C   <- matrix(0,nrow=ngroup,ncol=ngroup) # contact/mixing matrix of ngroup by ngroup
      
      C[1,1]  <- cr[1]*rho_unvac_unvac # contacts between unvac and unvac
      C[1,2]  <- cr[1]*rho_unvac_vac   # contacts between unvac and vac
      C[2,1]  <- cr[2]*rho_vac_unvac   # contacts between vac and unvac
      C[2,2]  <- cr[2]*rho_vac_vac*contact.increase   # contacts between vac and vac
      
      
      #SPECIFIC MODEL SET-UP (SET VE)
      
      sigma<-1-ve.infect #sigma is the vaccine efficacy for infectiousness 
      
      #Here we set alpha (vaccine efficacy against suscept.) depending on Vaccine Type
      if(vaccine.type=="AllOrNothing"){
        alpha=(1-0) 
        R_init<- c(0,Ntot[2]*ve.susc)  
        
      }else if(vaccine.type=="Leaky"){
        alpha=(1-ve.susc)
        R_init<- c(0,0)  
      }
      
      S_init    <- Ntot - I_init - R_init   # initial number susceptible, subtracting number infectious in each group
      CumFOI_init <- c(beta,beta*alpha)*C %*% as.matrix(I_init/Ntot) ####double check this
      FOI_init <- c(beta,beta*alpha)*C %*% as.matrix(I_init/Ntot)
      
      rho=c(beta,beta*alpha) 
      sigma.matrix=matrix(c(1,1,sigma,sigma),nrow=2,ncol=2)
      Reff_init<-as.matrix(S_init/Ntot)*rho*(as.matrix(C)*sigma.matrix)%*%as.matrix(c(1,1))*1/gamma
      
      
      CumIncid_init <- I_init
      Incid_init <- I_init
      P.time.infected_init<-rep(0,ngroup)
      
      parms <-c(gamma=gamma,
                beta=beta, 
                alpha=alpha,
                sigma=sigma,
                C=C)
      
      inits <-c(S = S_init,               
                I = I_init,           
                R = R_init,               
                CumIncid = CumIncid_init,
                Incid = Incid_init,
                CumFOI = CumFOI_init,
                FOI = FOI_init,
                P.time.infected=P.time.infected_init,
                Reff=Reff_init)  
      
      
      times = seq(0,300,1)  #run model for 300 days
      
      #RUN MODEL
      m.output= as.data.frame(lsoda(inits,times,ve_model_contact,parms))  
      
      #Calculating Cumulative Incidence OR and RR
      m.output$a<-m.output$CumIncid2
      m.output$b<-Ntot[2] -m.output$CumIncid2 
      m.output$c<-m.output$CumIncid1
      m.output$d<-Ntot[1] -m.output$CumIncid1
      
      m.output$OR.overall<-1-((m.output$a*m.output$d)/(m.output$b*m.output$c))
      m.output$RR.overall<-1-((m.output$a/(m.output$a + m.output$b))/(m.output$c/(m.output$c + m.output$d)))
      
      #Calculating Incidence OR and RR
      m.output$e<-m.output$Incid2
      m.output$f<-Ntot[2] -m.output$Incid2
      m.output$g<-m.output$Incid1
      m.output$h<-Ntot[1] -m.output$Incid1
      
      m.output$OR.perday<-1 - ((m.output$e*m.output$h)/(m.output$f*m.output$g))
      m.output$RR.perday<-1-((m.output$e/(m.output$e + m.output$f))/(m.output$g/(m.output$g + m.output$h)))
      
      #GET MINIMUM VALUES
      Contour.VE.df$min.OR.overall[nrow]<-min(m.output$OR.overall)
      Contour.VE.df$min.RR.overall[nrow]<-min(m.output$RR.overall)
      
      Contour.VE.df$min.OR.perday[nrow]<-min(m.output$OR.perday)
      Contour.VE.df$min.RR.perday[nrow]<-min(m.output$RR.perday)
      
      Contour.VE.df$ve.susc[nrow]<-ve.susc
      Contour.VE.df$ve.infect[nrow]<-ve.infect
      Contour.VE.df$contact.increase[nrow]<-contact.increase
      
      contact.increase<-contact.increase + 0.02
    }
    
    ve.susc<-ve.susc + 0.02
  }
  return(Contour.VE.df)
}

#Generating scenarios (can take a while to run)
VeS_Contact.contour.df<-veS.contact.contour(vaccine.type="AllOrNothing")

#write.csv(VeS_Contact.contour.df,"VES_and_ContactIncrease.csv")