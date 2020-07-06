# **Diagnosis RomeIII questionnaire**

#### Version 1: 21st July 2016 Ettje Tigchelaar
#### Version 2: 20th March 2020 Trishla Sinha 
#### Update on 25th March 2020 by Trishla Sinha: translating, editing scripts as per new requirements 

#### Based on the scoring document for the Rome foundation:  RomeIII_functional bowel disorders questionnaire.pdf
#### Execution of the code below will results in 29 additional variables:
   - 26 scoring variables representing the 26 questions of the questionnaire
   - variable IBS stating yes or no for IBS
   - variable Subtype specifying the IBS subtype C/D/M/U
   - variable FGID listing all functional gastrointestinal disorders from this questionnaire: IBS, IBS strict, functional constipation, functional diarrhea, functional bloating

#### Name of dataframe: RO 
If containing no sex information: Add a column sex to this data R0$sex and label it male and female
#### Variable names questions: ROME1 to ROME26 
#### score questions
#### new variables created: ROME1s - ROME26s
for(i in c("ROME1","ROME19","ROME21","ROME23","ROME25")){
  RO$j <- ifelse(RO[,i]=="never",0,
         ifelse(RO[,i]=="less than 1 day a month",1,
                ifelse(RO[,i]=="1 day a month",2,
                       ifelse(RO[,i]=="two to three days a month",3,
                              ifelse(RO[,i]=="1 day a week",4,
                                     ifelse(RO[,i]=="more than 1 day a week",5,
                                            ifelse(RO[,i]=="every day",6,NA)))))))
  colnames(RO)[colnames(RO)=='j'] <- c(paste(colnames(RO[i]),"s",sep=""))
}

RO$ROME2s <- ifelse(RO$ROME2=="no",0,
                 ifelse(RO$ROME2=="yes",1,
                        ifelse(RO$ROME2=="not applicable,because I do not menstruate (any more)",2,NA)))     

for(i in c("ROME3","ROME15","ROME17","ROME18","ROME20","ROME22","ROME24","ROME26")){
  RO$j <- ifelse(RO[,i]=="no",0,
                 ifelse(RO[,i]=="yes",1,NA))
  colnames(RO)[colnames(RO)=='j'] <- c(paste(colnames(RO[i]),"s",sep=""))
}

for(i in c("ROME4","ROME5","ROME6","ROME7","ROME8","ROME9","ROME10","ROME11","ROME12","ROME13","ROME14","ROME16")){
  RO$j <- ifelse(RO[,i]=="rarely or never",0,
                 ifelse(RO[,i]=="sometimes",1,
                        ifelse(RO[,i]=="often",2,
                               ifelse(RO[,i]=="most of the time",3,
                                      ifelse(RO[,i]=="always",4,NA)))))
  colnames(RO)[colnames(RO)=='j'] <- c(paste(colnames(RO[i]),"s",sep=""))
}


#### diagnostic criteria IBS
  a = ifelse(RO$ROME3s==1,1,0)                                                      # symptom onset 6 months prior to diagnosis
  b = ifelse(RO$ROME1s>2,1,0)                                                       # pain/discomfort (p/d) at least 2-3 days/month
  c = ifelse(RO$ROME2s==0 | RO$ROME2s==2 |(is.na(RO$ROME2s) & RO$sex=="male"),1,0) # pain not associated to menstrual bleeding
  d = ifelse(RO$ROME4s>0,1,0)                                                       # p/d better after bowel movement
  e = ifelse(RO$ROME5s>0 | RO$ROME6s>0,1,0)                                            # p/d associated with more or fewer stools
  f = ifelse(RO$ROME7s>0 | RO$ROME8s>0,1,0)                                            # p/d associated with looser/harder stools
  g = ifelse(RO$ROME1s>4,1,0)                                                       # stricter criteria for clinical trials: p/d >1 day/week
  h = ifelse(RO$ROME10s>0 & RO$ROME16s==0,1,0)                                         # criteria for IBS-C
  i = ifelse(RO$ROME10s==0 & RO$ROME16s>0,1,0)                                         # criteria for IBS-D
  j = ifelse(RO$ROME10s>0 & RO$ROME16s>0,1,0)                                          # criteria for IBS-M
  k = ifelse(RO$ROME10s==0 & RO$ROME16s==0,1,0)                                        # criteria for IBS-U

#### creating a variable for IBS yes/no and subtype of IBS
RO$IBS <- ifelse(a==1 & b==1 & c==1 & (d+e+f)>=2,"yes","no")
RO$Subtype <- ifelse(RO$IBS=="yes" & h==1,"C",
                     ifelse(RO$IBS=="yes" & i==1,"D",
                            ifelse(RO$IBS=="yes" & j==1,"M",
                                   ifelse(RO$IBS=="yes" & k==1,"U",RO$IBS))))

#### diagnostic criteria functional bloating/constipation/diarrhea

  l = ifelse(RO$ROME20s==1 & RO$ROME19s>2 & (RO$ROME13s<5 | RO$ROME14s==0) & (RO$ROME15s<5 | RO$ROME16s==0) & (RO$ROME17s<5 | RO$ROME18s==0),1,0) # functional bloating
  m = ifelse(RO$ROME15s==1 & RO$ROME7s==0,1,0)            # functional constipation
  m1 = ifelse(RO$ROME11s>1,1,0)
  m2 = ifelse(RO$ROME10s>1,1,0)
  m3 = ifelse(RO$ROME12s>0,1,0)
  m4 = ifelse(RO$ROME13s>0,1,0)
  m5 = ifelse(RO$ROME14s>0,1,0)
  m6 = ifelse(RO$ROME9s>1,1,0)
  n = ifelse(RO$ROME18s==1 & RO$ROME17s==1 & RO$ROME1s==0,1,0)  # functional diarrhea

#### creating a variable for all functional gastrointestinal disorders including IBS strict, IBS, functional constipation, functional bloating and functional diarrhea
RO$FGID <- ifelse(RO$IBS=="yes" & g==1, "IBS strict",
                  ifelse(RO$IBS=="yes","IBS",
                      ifelse(RO$IBS=="no" & m==1 & (m1+m2+m3+m4+m5+m6)>=2,"functional constipation",
                            ifelse(RO$IBS=="no" & n==1,"functional diarrhea",
                                   ifelse(RO$IBS=="no" & l==1,"functional bloating",RO$IBS)))))


##### Note that IBS strict criteria is not used for this pupose as its a stricter definition for clinical trials 
##### In FGID, NA's to No's as all participants answere the questionaires 
