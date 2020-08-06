
#Loading package and data
library(survival)
clinical.data<-data.frame(read.table("BRCA_clinical.txt",header=TRUE,stringsAsFactors=FALSE))
Mutation.data<-data.frame(read.table("Annotated_Somatic_mutation.txt",header=TRUE,stringsAsFactors=FALSE,fill=TRUE))


#Selecting TP53 Mutated Sample Barcode 
CaseBarcode<-unique(Mutation.data[which(Mutation.data$Hugo_Symbol=="TP53"&!(Mutation.data$Variant_Classification %in% c("In_Frame_Del","Silent"))),])$Sample_Barcode
#CaseBarcode<-unique(Mutation.data[which(Mutation.data$Hugo_Symbol=="TP53"),])$Sample_Barcode


#Dividing Cohort into Case and Control group
Case.Followup<-data.frame(clinical.data[which(clinical.data$bcr_patient_barcode %in% CaseBarcode),c(2,14,15,16)],group="Case",stringsAsFactors=FALSE)
Control.Followup<-data.frame(clinical.data[which(!(clinical.data$bcr_patient_barcode %in% CaseBarcode)),c(2,14,15,16)],group="Control",stringsAsFactors=FALSE)
str(Control.Followup)


#Unifying data
Followup.data<-rbind.data.frame(Case.Followup,Control.Followup)
Followup.data[which(Followup.data$vital_status=="Dead"),3]<-Followup.data[which(Followup.data$vital_status=="Dead"),4]
str(Followup.data)

#Converting character data into numeric data
Followup.data$days_to_last_followup<-as.numeric(Followup.data$days_to_last_followup)

#Survival Analysis
logtest.byTP53 <- survdiff(Surv(days_to_last_followup, vital_status=="Dead")~as.factor(group),data=Followup.data)
log.p<-1 - pchisq(logtest.byTP53$chisq, 1)
log.p
surv.byTP53 <- survfit(Surv(days_to_last_followup, vital_status=="Dead")~as.factor(group),data=Followup.data)

#Drawing Kaplan-Meier Survival Plot
plot(surv.byTP53,col=c("red","black"),main=paste("Survival Analysis - TP53(p=",round(log.p,digit=4),")"),xlab="Days of Follow up",
     ylab="Overall Survival(%)",mark.time=TRUE)


#Adding Legend

legend("bottomleft",legend=c(paste("Case:",254,"/",289),
                             paste("Control:",728,"/",796)),
       col=c("red","black"),lty=c(1,1))

