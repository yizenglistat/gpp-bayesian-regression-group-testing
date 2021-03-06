# To remind 
# Z = A matrix of testing responses whose jth row is of the form (col1=Z_j, col2=cj, col3=assay used, col4:col(4+cj-1)=indices of the individuals assigned to the jth pool) 
# X = Covariate matrix in which the ith row contains the covariate information pertaining to the ith individual
# Y = matrix whose ith row is of the form (col1=0, col2=number of pools ith individual involved in, col3:col(3+col2-1)=pools ith individual was assigned to)
##################################################################################################################
# Functions that implement Individual testing, Master pool testing, Dorfman testing, Array testing


#############################################################################
# R function: This function simulates individual level testin and stores the 
#             testing responses in accordance to the data structure 
#             required to fit the Bayesian GT regression model presented
#             in McMahan et al. (2014+) via the R function Bayes.reg.GT
#
# Input: 
# Y_true = The true statuses of the individuals
# Se = The testing sensitivity used for both pools and individual testing
# Sp = The testing specificity used for both pools and individual testing
# 

# indTesting<-function(Y_true,Se,Sp){
# N<-length(Y_true)
# Y<-rbinom(N,1,(Se*Y_true+(1-Sp)*(1-Y_true)))
# return(list("Z"=Y, "Y"=Y))
# }

indTesting<-function(Y_true,Se,Sp,cj=1){
   N<-length(Y_true)
   Jmax<-N/cj 
   J<-1

   Y<-matrix(-99,nrow=N, ncol=3) 
   Z<-matrix(-99,nrow=Jmax,ncol=cj+3) 


   for(j in 1:(N/cj)){
   prob<-ifelse(sum(Y.true[((j-1)*cj+1):(j*cj)])>0,Se,1-Sp)
   Z[J,1]<-rbinom(n=1,size=1,prob=prob)
   Z[J,2]<-cj
   Z[J,3]<-1
   Z[J,4:(cj+3)]<-((j-1)*cj+1):(j*cj)
   Y[((j-1)*cj+1):(j*cj),2]<-1
   Y[((j-1)*cj+1):(j*cj),3]<-J
   J<-J+1
   }

   J<-J-1
   Z<-Z[1:J,]

   return(list("Z"=Z, "Y"=Y))
}

#####################################################################
# R function: This function simulates Initial pool testing and stores the 
#             testing responses in accordance to the data structure 
#             required to fit the Bayesian GT regression model presented
#             in McMahan et al. (2012+) via the R function Bayes.reg.GT
#
# Input: 
# Y_true = The true statuses of the individuals
# Se = A vector of testing sensitivities, where the first element is the
#      testing sensitivity for the pools and the second entry is the 
#      test sensitivity for individual testing
# Sp = A vector of testing specificities, where the first element is the
#      testing specificity for the pools and the second entry is the 
#      test specificity for individual testing
# cj = pool size to be used (Note: The number of individuals N should be 
#      evenly divisible by cj, this is only for decoding purposes; i.e., 
#      the regression methods do not require this condition)


masterpoolTesting<-function(Y_true,Se,Sp,cj){
   N<-length(Y_true)
   Jmax<-N/cj 
   J<-1

   Y<-matrix(-99,nrow=N, ncol=3) 
   Z<-matrix(-99,nrow=Jmax,ncol=cj+3) 


   for(j in 1:(N/cj)){
      prob<-ifelse(sum(Y_true[((j-1)*cj+1):(j*cj)])>0,Se,1-Sp)
      Z[J,1]<-rbinom(n=1,size=1,prob=prob)
      Z[J,2]<-cj
      Z[J,3]<-1
      Z[J,4:(cj+3)]<-((j-1)*cj+1):(j*cj)
      Y[((j-1)*cj+1):(j*cj),2]<-1
      Y[((j-1)*cj+1):(j*cj),3]<-J
      J<-J+1
   }

   J<-J-1
   Z<-Z[1:J,]

   return(list("Z"=Z, "Y"=Y))
}



#####################################################################
# R function: This function simulates Dorfman decoding and stores the 
#             testing responses in accordance to the data structure 
#             required to fit the Bayesian GT regression model presented
#             in McMahan et al. (2012+) via the R function Bayes.reg.GT
#
# Input: 
# Y_true = The true statuses of the individuals
# Se = A vector of testing sensitivities, where the first element is the
#      testing sensitivity for the pools and the second entry is the 
#      test sensitivity for individual testing
# Sp = A vector of testing specificities, where the first element is the
#      testing specificity for the pools and the second entry is the 
#      test specificity for individual testing
# cj = pool size to be used (Note: The number of individuals N should be 
#      evenly divisible by cj, this is only for decoding purposes; i.e., 
#      the regression methods do not require this condition)


dorfmanTesting_diff_cj<-function(Y_true,Se,Sp,cj){
   N<-length(Y_true)
   Jmax<-N+N/cj
   J<-1

   Y<-matrix(-99,nrow=N, ncol=4)
   Z<-matrix(-99,nrow=Jmax,ncol=cj+3)


   for(j in 1:length(cj)){

      prob<-ifelse(sum(Y_true[((j-1)*cj[j]+1):(j*cj[j])])>0,Se[[as.character(cj[j])]],1-Sp[[as.character(cj[j])]])
      Z[J,1]<-rbinom(n=1,size=1,prob=prob)
      Z[J,2]<-cj[j]
      Z[J,3]<-1   ## swab pool used
      Z[J,4:(cj[j]+3)]<-((j-1)*cj[j]+1):(j*cj[j])
      Y[((j-1)*cj[j]+1):(j*cj[j]),2]<-1
      Y[((j-1)*cj[j]+1):(j*cj[j]),3]<-J
      J<-J+1
      if(Z[J-1,1]==1){
         for(k in ((j-1)*cj[j]+1):(j*cj[j])){
            prob<-ifelse(Y_true[k]>0,Se[[as.character(1)]],1-Sp[[as.character(1)]])
            Z[J,1]<- rbinom(n=1,size=1,prob=prob)
            Z[J,2]<-1
            Z[J,3]<-2   ## swab ind used
            Z[J,4]<-k
            Y[k,2]<-2
            Y[k,4]<-J
            J<-J+1
         }
      }
   }

   J<-J-1
   Z<-Z[1:J,]

   return(list("Z"=Z, "Y"=Y))
}

dorfmanTesting_diff<-function(Y_true,Se,Sp,cj){
   N<-length(Y_true)
   Jmax<-N+N/cj
   J<-1

   Y<-matrix(-99,nrow=N, ncol=4)
   Z<-matrix(-99,nrow=Jmax,ncol=cj+3)


   for(j in 1:(N/cj)){
      prob<-ifelse(sum(Y_true[((j-1)*cj+1):(j*cj)])>0,Se[1],1-Sp[1])
      Z[J,1]<-rbinom(n=1,size=1,prob=prob)
      Z[J,2]<-cj
      Z[J,3]<-1   ## swab pool used
      Z[J,4:(cj+3)]<-((j-1)*cj+1):(j*cj)
      Y[((j-1)*cj+1):(j*cj),2]<-1
      Y[((j-1)*cj+1):(j*cj),3]<-J
      J<-J+1
      if(Z[J-1,1]==1){
         for(k in ((j-1)*cj+1):(j*cj)){
            prob<-ifelse(Y_true[k]>0,Se[2],1-Sp[2])
            Z[J,1]<- rbinom(n=1,size=1,prob=prob)
            Z[J,2]<-1
            Z[J,3]<-2   ## swab ind used
            Z[J,4]<-k
            Y[k,2]<-2
            Y[k,4]<-J
            J<-J+1
         }
      }
   }

   J<-J-1
   Z<-Z[1:J,]

   return(list("Z"=Z, "Y"=Y))
}


#####################################################################
# R function: This function simulates Dorfman decoding and stores the 
#             testing responses in accordance to the data structure 
#             required to fit the Bayesian GT regression model presented
#             in McMahan et al. (2012+) via the R function Bayes.reg.GT
#
# Input: 
# Y_true = The true statuses of the individuals
# Se = The testing sensitivity used for both pools and individual testing
# Sp = The testing specificity used for both pools and individual testing
# cj = pool size to be used (Note: The number of individuals N should be 
#      evenly divisible by cj, this is only for decoding purposes; i.e., 
#      the regression methods do not require this condition)



dorfmanTesting_same<-function(Y_true,Se,Sp,cj){
N<-length(Y_true)
Jmax<-N+N/cj 
J<-1

Y<-matrix(-99,nrow=N, ncol=4) 
Z<-matrix(-99,nrow=Jmax,ncol=cj+3) 


for(j in 1:(N/cj)){
prob<-ifelse(sum(Y_true[((j-1)*cj+1):(j*cj)])>0,Se,1-Sp)
Z[J,1]<-rbinom(n=1,size=1,prob=prob)
Z[J,2]<-cj
Z[J,3]<-1
Z[J,4:(cj+3)]<-((j-1)*cj+1):(j*cj)
Y[((j-1)*cj+1):(j*cj),2]<-1
Y[((j-1)*cj+1):(j*cj),3]<-J
J<-J+1
if(Z[J-1,1]==1){
for(k in ((j-1)*cj+1):(j*cj)){
prob<-ifelse(Y_true[k]>0,Se,1-Sp)
Z[J,1]<- rbinom(n=1,size=1,prob=prob)
Z[J,2]<-1
Z[J,3]<-1
Z[J,4]<-k
Y[k,2]<-2
Y[k,4]<-J
J<-J+1
}
}
}

J<-J-1
Z<-Z[1:J,]

return(list("Z"=Z, "Y"=Y))
}


#####################################################################
# R function: This function simulates Array decoding and stores the 
#             testing responses in accordance to the data structure 
#             required to fit the Bayesian GT regression model presented
#             in McMahan et al. (2012+) via the R function Bayes.reg.GT
#
# Input: 
# Y_true = The true statuses of the individuals
# Se = A vector of testing sensitivities, where the first element is the
#      testing sensitivity for the pools and the second entry is the 
#      test sensitivity for individual testing
# Sp = A vector of testing specificities, where the first element is the
#      testing specificity for the pools and the second entry is the 
#      test specificity for individual testing
# cj = row and column pool sizes to be used (Note: The number of individuals 
#      N should be evenly divisible by cj^2, this is only for decoding 
#      purposes; i.e., the regression methods do not require this condition)


arrayTesting_diff<-function(Y_true, Se, Sp, cj){

N<-length(Y_true)
Jmax<-2*N/cj +N
J<-1
AT<-N/(cj^2)

Y<-matrix(-99,nrow=N, ncol=5) 
Z<-matrix(-99,nrow=Jmax,ncol=cj+3) 

Y.A<-array(-99,c(cj,cj,AT))
ID.A<-array(-99,c(cj,cj,AT))
ind<-1
for(i in 1:AT){
for(j in 1:cj){
for(m in 1:cj){
Y.A[m,j,i]<-Y_true[ind]
ID.A[m,j,i]<-ind
ind<-ind+1
}}}



for(s in 1:AT){

array.yk<-Y.A[,,s]
array.id<-ID.A[,,s]

a<-rep(0,nrow(array.yk))
b<-rep(0,ncol(array.yk))

for(i in 1:cj){
   prob<- ifelse(sum(array.yk[i,])>0, Se[1], 1-Sp[1])
   g<- rbinom(n=1,size=1,prob=prob)
   a[i]<-g
   Z[J,1]<-g 
   Z[J,2]<-cj 
   Z[J,3]<-1
   Z[J,4:(cj+3)]<-array.id[i,]
   Y[array.id[i,],2]<-2
   Y[array.id[i,],3]<-J
   J<-J+1
}
for(j in 1:cj){
   prob<- ifelse(sum(array.yk[,j])>0, Se[1], 1-Sp[1])
   g<- rbinom(n=1,size=1,prob=prob)
   b[j]<-g
   Z[J,1]<-g 
   Z[J,2]<-cj 
   Z[J,3]<-1
   Z[J,4:(cj+3)]<-array.id[,j]
   Y[array.id[,j],4]<-J
   J<-J+1
}


if(sum(a)>0 & sum(b)>0){
array.yk1<-as.matrix(array.yk[(a==1),(b==1)])
array.id1<-as.matrix(array.id[(a==1),(b==1)])
for(i in 1:nrow(array.yk1)){
for(j in 1:ncol(array.yk1)){
   prob<- ifelse(array.yk1[i,j]>0, Se[2], 1-Sp[2])
   g<- rbinom(n=1,size=1,prob=prob)
   Z[J,1]<-g 
   Z[J,2]<-1 
   Z[J,3]<-2
   Z[J,4]<-array.id1[i,j]
   Y[array.id1[i,j],2]<-3
   Y[array.id1[i,j],5]<-J
   J<-J+1
}}}



if(sum(a)>0 & sum(b)==0){
array.yk1<-as.matrix(array.yk[(a==1),])
array.id1<-as.matrix(array.id[(a==1),])
for(i in 1:nrow(array.yk1)){
for(j in 1:ncol(array.yk1)){
   prob<- ifelse(array.yk1[i,j]>0, Se[2], 1-Sp[2])
   g<- rbinom(n=1,size=1,prob=prob)
   Z[J,1]<-g 
   Z[J,2]<-1 
   Z[J,3]<-2
   Z[J,4]<-array.id1[i,j]
   Y[array.id1[i,j],2]<-3
   Y[array.id1[i,j],5]<-J
   J<-J+1
}}}

if(sum(a)==0 & sum(b)>0){
array.yk1<-as.matrix(array.yk[,(b==1)])
array.id1<-as.matrix(array.id[,(b==1)])
for(i in 1:nrow(array.yk1)){
for(j in 1:ncol(array.yk1)){
   prob<- ifelse(array.yk1[i,j]>0, Se[2], 1-Sp[2])
   g<- rbinom(n=1,size=1,prob=prob)
   Z[J,1]<-g 
   Z[J,2]<-1 
   Z[J,3]<-2
   Z[J,4]<-array.id1[i,j]
   Y[array.id1[i,j],2]<-3
   Y[array.id1[i,j],5]<-J
   J<-J+1
}}}

} 

J<-J-1
Z<-Z[1:J,]

return(list("Z"=Z, "Y"=Y))
} #End function




#####################################################################
# R function: This function simulates Array decoding and stores the 
#             testing responses in accordance to the data structure 
#             required to fit the Bayesian GT regression model presented
#             in McMahan et al. (2012+) via the R function Bayes.reg.GT
#
# Input: 
# Y_true = The true statuses of the individuals
# Se = The testing sensitivity used for both pools and individual testing
# Sp = The testing specificity used for both pools and individual testing
# cj = row and column pool sizes to be used (Note: The number of individuals 
#      N should be evenly divisible by cj^2, this is only for decoding 
#      purposes; i.e., the regression methods do not require this condition)


arrayTesting_same<-function(Y_true, Se, Sp, cj){

N<-length(Y_true)
Jmax<-2*N/cj +N
J<-1
AT<-N/(cj^2)

Y<-matrix(-99,nrow=N, ncol=5) 
Z<-matrix(-99,nrow=Jmax,ncol=cj+3) 

Y.A<-array(-99,c(cj,cj,AT))
ID.A<-array(-99,c(cj,cj,AT))
ind<-1
for(i in 1:AT){
for(j in 1:cj){
for(m in 1:cj){
Y.A[m,j,i]<-Y_true[ind]
ID.A[m,j,i]<-ind
ind<-ind+1
}}}



for(s in 1:AT){

array.yk<-Y.A[,,s]
array.id<-ID.A[,,s]

a<-rep(0,nrow(array.yk))
b<-rep(0,ncol(array.yk))

for(i in 1:cj){
   prob<- ifelse(sum(array.yk[i,])>0, Se, 1-Sp)
   g<- rbinom(n=1,size=1,prob=prob)
   a[i]<-g
   Z[J,1]<-g 
   Z[J,2]<-cj 
   Z[J,3]<-1
   Z[J,4:(cj+3)]<-array.id[i,]
   Y[array.id[i,],2]<-2
   Y[array.id[i,],3]<-J
   J<-J+1
}
for(j in 1:cj){
   prob<- ifelse(sum(array.yk[,j])>0, Se, 1-Sp)
   g<- rbinom(n=1,size=1,prob=prob)
   b[j]<-g
   Z[J,1]<-g 
   Z[J,2]<-cj 
   Z[J,3]<-1
   Z[J,4:(cj+3)]<-array.id[,j]
   Y[array.id[,j],4]<-J
   J<-J+1
}


if(sum(a)>0 & sum(b)>0){
array.yk1<-as.matrix(array.yk[(a==1),(b==1)])
array.id1<-as.matrix(array.id[(a==1),(b==1)])
for(i in 1:nrow(array.yk1)){
for(j in 1:ncol(array.yk1)){
   prob<- ifelse(array.yk1[i,j]>0, Se, 1-Sp)
   g<- rbinom(n=1,size=1,prob=prob)
   Z[J,1]<-g 
   Z[J,2]<-1 
   Z[J,3]<-1
   Z[J,4]<-array.id1[i,j]
   Y[array.id1[i,j],2]<-3
   Y[array.id1[i,j],5]<-J
   J<-J+1
}}}



if(sum(a)>0 & sum(b)==0){
array.yk1<-as.matrix(array.yk[(a==1),])
array.id1<-as.matrix(array.id[(a==1),])
for(i in 1:nrow(array.yk1)){
for(j in 1:ncol(array.yk1)){
   prob<- ifelse(array.yk1[i,j]>0, Se, 1-Sp)
   g<- rbinom(n=1,size=1,prob=prob)
   Z[J,1]<-g 
   Z[J,2]<-1 
   Z[J,3]<-1
   Z[J,4]<-array.id1[i,j]
   Y[array.id1[i,j],2]<-3
   Y[array.id1[i,j],5]<-J
   J<-J+1
}}}

if(sum(a)==0 & sum(b)>0){
array.yk1<-as.matrix(array.yk[,(b==1)])
array.id1<-as.matrix(array.id[,(b==1)])
for(i in 1:nrow(array.yk1)){
for(j in 1:ncol(array.yk1)){
   prob<- ifelse(array.yk1[i,j]>0, Se, 1-Sp)
   g<- rbinom(n=1,size=1,prob=prob)
   Z[J,1]<-g 
   Z[J,2]<-1 
   Z[J,3]<-1
   Z[J,4]<-array.id1[i,j]
   Y[array.id1[i,j],2]<-3
   Y[array.id1[i,j],5]<-J
   J<-J+1
}}}

} 

J<-J-1
Z<-Z[1:J,]

return(list("Z"=Z, "Y"=Y))
} #End function





#####################################################################
# R function: This function translates the real data of Chlamydia test  
#             to the format for model fitting. Output will be matrix
#             X, Y and Z
#             
#
# Input: 
# Real Data

format_realdata <- function(data.new){

   pool.id  <- as.numeric(data.new$Pool.ID)
   data.gts <- data.new[pool.id!="NaN",]
   data.ind <- data.new[pool.id=="NaN",]

   ##################################################################
   # Creating the design matrices
   X.ind<-cbind(
                      as.numeric(as.character(data.ind$Age)),  
                      as.numeric(data.ind$Race=="W"), 
                      as.numeric(data.ind$Risk.New.Partner=="Y"),
                      as.numeric(data.ind$Risk.Multiple.Partner=="Y"),
                      as.numeric(data.ind$Risk.Contact=="Y"),
                      as.numeric(data.ind$Symptom=="Y"),
                      as.numeric(data.ind$Specimen.Type=="Swab")
   )


   X.pool<-cbind(
                      as.numeric(as.character(data.gts$Age)),  
                      as.numeric(data.gts$Race=="W"), 
                      as.numeric(data.gts$Risk.New.Partner=="Y"),
                      as.numeric(data.gts$Risk.Multiple.Partner=="Y"),
                      as.numeric(data.gts$Risk.Contact=="Y"),
                      as.numeric(data.gts$Symptom=="Y"),
                      as.numeric(data.gts$Specimen.Type=="Swab")
   )

   dim(X.ind)
   X.ind[1:10,]
   ###################################################################
   # Building the Z and Y matrice

   Z.ind.C<-matrix(-99,nrow=dim(data.ind)[1],ncol=8)
   Y.ind.C<-matrix(-99,nrow=dim(data.ind)[1],ncol=4)

   Z.pool.C<-matrix(-99,nrow=dim(data.gts)[1],ncol=8)
   Y.pool.C<-matrix(-99,nrow=dim(data.gts)[1],ncol=4)


   ##################################
   # For GT testing

   id.ind<-1:(dim(data.gts)[1])

   pool.id<-unique(data.gts$Pool.ID)
   n<-length(pool.id)

   track.CT<-1


   #############################
   #############################
   ### For loop starts

   for(i in 1:n){

   temp<-data.gts[data.gts$Pool.ID==pool.id[i],]
   temp.id<-id.ind[data.gts$Pool.ID==pool.id[i]]

   CT.res <- as.numeric(temp$CT.Result=="P")

   cj<-length(CT.res)

   CT.retest<- sum(CT.res)>0

   Z.pool.C[track.CT,1]<-CT.retest
   Z.pool.C[track.CT,2]<-cj
   Z.pool.C[track.CT,3]<-1 # Swab Assay
   Z.pool.C[track.CT,4:(cj+3)]<-temp.id

   Y.pool.C[temp.id,1]<-0
   Y.pool.C[temp.id,2]<-1
   Y.pool.C[temp.id,3]<-track.CT


   if(CT.retest==0){
   track.CT<-track.CT+1
   }

   if(CT.retest>0){
   tid<-(track.CT+1):(track.CT+cj)
   Z.pool.C[tid,1]<-CT.res
   Z.pool.C[tid,2]<-1
   Z.pool.C[tid,3]<-1 # Swab Assay
   Z.pool.C[tid,4]<-temp.id

   Y.pool.C[temp.id,1]<-CT.res
   Y.pool.C[temp.id,2]<-2
   Y.pool.C[temp.id,4]<-tid
   track.CT<-max(tid)+1
   }


   }

   ### For loop ends
   #############################
   #############################


   Z.pool.C<-Z.pool.C[1:(track.CT-1),]

   ZnC<-dim(Z.pool.C)[1]
   YnC<-dim(Y.pool.C)[1]


   #################################
   # For individual level testing 

   Z.ind.C[,1]<-as.numeric(data.ind$CT.Result=="P") # Test outcome for Chlamydia
   Z.ind.C[,2]<-1                                   # Pool size    
   Z.ind.C[,3]<- as.numeric(data.ind$Specimen.Type=="Urine")+1 #Urine Assay
   Z.ind.C[,4]<-(1:(dim(data.ind)[1]))+YnC


   Y.ind.C[,1]<-as.numeric(data.ind$CT.Result=="P") # Initial true status
   Y.ind.C[,2]<-1
   Y.ind.C[,3]<-1:(dim(data.ind)[1])+ZnC


   ###################################################
   # Putting everything together 

   X<-rbind(X.pool,X.ind)
   colnames(X) <- c("Age","Race_W","Risk.New.Partner","Risk.Multiple.Partner","Risk.Contact",
   "Symptom","Specimen.Type")

   Z.C<-rbind(Z.pool.C,Z.ind.C)
   Y.C<-rbind(Y.pool.C,Y.ind.C)

  return(list(X=X,Y=Y.C,Z=Z.C))

}























