### survival function
multicoxfull<-function(data, entrytime, score, outcome){
       if (!outcome %in% c('OS', 'DOC', 'CSS')){
         stop("Invalid outcome. Please choose 'OS', 'DOC', or 'CSS'.")
       }
  
       fml_base <- paste("Surv(time = ", entrytime, ", time2 = timeD, ",
                            if (outcome == 'OS') "death_all)" else
                              if (outcome =='DOC') "death_crccp==2)" else "death_crccp==1)", 
                               "~Age_diag+Sex+TNM_adj+BMI+high_blood_pressure+NSAIDS+smoking+cardiovad+statins+chemradther+crc2sites+physical_activity_lifeav+evercol+hormone_replace+ ")
                   
       fml <- formula(paste(fml_base, score))               
         
       model <- coxph(fml,data = data)
         
       b<-as.data.frame(ShowRegTable(model, exp = TRUE, digits = 2, pDigits = 3,
                                     printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T)
       
       
       score_rows <- grep(score, rownames(b))
       
       b<-b[score_rows, ]
       colnames(b)[c(1, 2)]<-c('aHR(95% CI)','aPvalue')
       b$score<-rownames(b)
       b$outcome<-outcome
       return(b)
}


multicoxpart<-function(data, entrytime, score, outcome){
  if (!outcome %in% c('OS', 'DOC', 'CSS')){
    stop("Invalid outcome. Please choose 'OS', 'DOC', or 'CSS'.")
  }
  
  fml_base <- paste("Surv(time = ", entrytime, ", time2 = timeD, ",
                    if (outcome == 'OS') "death_all)" else
                      if (outcome =='DOC') "death_crccp==2)" else "death_crccp==1)", 
                     #"~Age_diag+Sex+TNM_adj+BMI+high_blood_pressure+NSAIDS+smoking+cardiovad+statins+chemradther+crc2sites+physical_activity_lifeav+evercol+hormone_replace+ ")
                     "~Age_diag+Sex+TNM_adj+ ")

  fml <- formula(paste(fml_base, score))               
  
  model <- coxph(fml,data = data)
  
  b<-as.data.frame(ShowRegTable(model, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T)
  
  
  score_rows <- grep(score, rownames(b))
  
  b<-b[score_rows, ]
  colnames(b)[c(1, 2)]<-c('aHR(95% CI)','aPvalue')
  b$score<-rownames(b)
  b$outcome<-outcome
  return(b)
}





