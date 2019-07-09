### This code corresponds to the review article "Life expectancy -
### frequently used, but hardly understood"

### Clears workspace
rm(list=ls())

### This package is required to run the code
library(dplyr)

##############################
######### Functions ##########
##############################

TMR.MAD.CAL.function <- function(HMD.frame, Period) {

    qx.frame <- HMD.frame[,c("Year", "Age", "qx")]
    qx.frame$Cohort <- qx.frame$Year-qx.frame$Age
    qx.frame$px <- 1-qx.frame$qx
    max.age <- max(qx.frame$Age)
    min.age <- 30 # Measures refer to ages 30+
    data.select.cohorts <- filter(qx.frame, Cohort>=(Period-max.age) &
                                            Cohort<=Period)
    data.select.age <- filter(data.select.cohorts, Age>=min.age)
    data.select <- arrange(data.select.age, Cohort)
    data.select <- data.select %>%
                           group_by(Cohort)  %>%
    mutate(lx=c(1,cumprod(px)[-length(px)])) %>%
    mutate(pxC=cumprod(px))

    TMR.lx <- c()

    for (i in 1:81) {
            i.age <- i+29
            TMR.lx[i] <- data.select$lx[data.select$Cohort==(Period-i.age)& data.select$Age==(i.age)]
            }

    TMR.qx <- qx.frame$qx[qx.frame$Year==Period & qx.frame$Age>=30]

    TMR <- sum(TMR.lx*TMR.qx)

    TMR.dx <- TMR.qx*TMR.lx
    age <- 30:110
    MAD <- (sum((age+0.5)*TMR.dx)/sum(TMR.dx))-30

    CALpx <- c()

    for (i in 1:81) {
     i.age <- i+29
     CALpx[i] <- data.select$pxC[data.select$Cohort==(Period-i.age) & data.select$Age==(i.age)]
      }

    CAL <- sum(CALpx)+0.5

    output <- data.frame(TMR=TMR,
                     MAD=MAD,
                     CAL=CAL
                     )
    return(output)
}

### Thanks to Marcus Ebeling, who provided this useful function during a
### mortality analysis class taught at the University of Rostock

life.table <- function(mx){
    ax <- c(0.14, rep(0.5, length(mx)-1))
    qx <- mx/(1+(1-ax)*mx)
    qx[length(qx)] <- 1
    qx[qx > 1] <- 1
    px <- 1-qx
    lx <- c(100000, (cumprod(px)*100000)[1:(length(px)-1)])
    dx <- c(-diff(lx), lx[length(lx)])
    Lx1 <- lx[-1]+ax[-length(ax)]*dx[-length(dx)]
    open.Lx <-  ifelse( mx[length(mx)] == 0, 0, dx[length(dx)]/mx[length(mx)])
    Lx <- c(Lx1, open.Lx)
    Tx <- rev(cumsum(rev(Lx)))
    ex <- Tx/lx

    return(data.frame(qx=qx, px = px, ax = ax, lx = lx , dx = dx, Lx= Lx,
                  Tx = Tx, ex = ex))
}


Tempo.adj.LE.function <- function(HMD.frame, Output.frame) {

    Tempo.LE <- c()
    LE30 <- c()

    for (i in 1:7) {
        i.year <- 2008+i
        select.year <- filter(HMD.frame, Year==i.year)
        LE30[i] <- select.year$ex[31]
        the.mx <- select.year$mx
        the.TMR <- Output.frame[,i]$TMR
        adj.mx.f <- c(the.mx[c(1:30)], the.mx[c(31:111)]/the.TMR)
        Tempo.LE[i] <- life.table(adj.mx.f)$ex[31]
        }

    add <- rbind(LE30, Tempo.LE)
    Output.frame <- rbind(Output.frame, add)
    rownames(Output.frame) <- c("TMR","MAD","CAL","PLE","LE*")
    colnames(Output.frame) <- paste("X",2009:2015, sep="")

    return(Output.frame)
}


#############################
####### Calculations ########
#############################


### All the calculations are based on HMD period life table data (www.mortality.org).
### Please note, results might deviate from the figures in the paper
### due to changes in the HMD database
### Our estimates are based on the files downloaded on the 4th of
### February 2019



############
###Belgium##
############

### folder with HMD data for Belgium (life table data for females)
setwd("s:/LETHE/Geron-Tempo/Belgium/")

HMD.Belgium <- read.table("fltper_1x1.txt", header=TRUE, skip=2)
HMD.Belgium$Age <- as.numeric(as.character(HMD.Belgium$Age))
HMD.Belgium$Age[is.na(HMD.Belgium$Age)] <- 110 # "110+" becomes 110
HMD.Belgium$qx <- as.numeric(as.character(HMD.Belgium$qx))
HMD.Belgium$mx <- as.numeric(as.character(HMD.Belgium$mx))
HMD.Belgium$ex <- as.numeric(as.character(HMD.Belgium$ex))

### Do not worry about the warning message. Belgium has a gap in the
### data between 1914-1918. However, this is not relevant for our estimations.

Belgium.output <- sapply(c(2009:2015), TMR.MAD.CAL.function, HMD.frame=HMD.Belgium)
Belgium <- Tempo.adj.LE.function(HMD.frame=HMD.Belgium, Output.frame=Belgium.output)

### Output table
Belgium


################
### Figure 1 ###
################

my.years <- 2009:2015

par(mar = c(5, 4, 4, 4) + 0.1)
plot(my.years, Belgium["PLE",], ylim=c(52.5,56.5),
     type="l", ylab="PLE", xlab="Calendar Year", main="Belgium")
par(new = TRUE)
plot(my.years, Belgium["TMR",], ylim=c(0.75,1), type
     = "l", axes = FALSE, bty = "n", xlab = "", ylab = "", lty=2)
axis(side=4, at = seq(0.8,1,0.05))
mtext("TMR", side=4, line=3)
legend("topleft", legend=c("PLE","TMR"), lty=c(1,2), bty="n")

################
### Figure 2 ###
################

change.LE30 <- c(as.numeric(Belgium["PLE",])-as.numeric(Belgium["PLE",])[1])
change.MAD <- c(as.numeric(Belgium["MAD",])-as.numeric(Belgium["MAD",])[1])
change.CAL <- c(as.numeric(Belgium["CAL",])-as.numeric(Belgium["CAL",])[1])
change.TempoLE <- c(as.numeric(Belgium["LE*",])-as.numeric(Belgium["LE*",])[1])

plot(my.years, change.LE30, ylim=c(0,1.5), type="l",
     xlab="Calendar Year", ylab="Change since 2009", lwd=2, main="Belgium")
lines(my.years, change.MAD, col="gray50", lty=2, lwd=2)
lines(my.years, change.CAL, col="gray50", lty=3, lwd=2)
lines(my.years, change.TempoLE, col="gray50", lty=1, lwd=2)
legend(2010, 1.2, legend=c("PLE","MAD","CAL","LE*"), lty=c(1,2,3,1),
       col=c("black",rep("gray50",3)), bty="n", ncol=2)


############
###France###
############

### folder with data for France
setwd("s:/LETHE/Geron-Tempo/France/")
HMD.France <- read.table("fltper_1x1.txt", header=TRUE, skip=2)
HMD.France$Age <- as.numeric(as.character(HMD.France$Age))
HMD.France$Age[is.na(HMD.France$Age)] <- 110


France.output <- sapply(c(2009:2015), TMR.MAD.CAL.function, HMD.frame=HMD.France)
France <- Tempo.adj.LE.function(HMD.frame=HMD.France, Output.frame=France.output)
France

################
### Figure 1 ###
################

my.years <- 2009:2015

par(mar = c(5, 4, 4, 4) + 0.1)
plot(my.years, France["PLE",], ylim=c(52.5,56.5),
     type="l", ylab="PLE", xlab="Calendar Year", main="France")
par(new = TRUE)
plot(my.years, France["TMR",], ylim=c(0.75,1), type
     = "l", axes = FALSE, bty = "n", xlab = "", ylab = "", lty=2)
axis(side=4, at = seq(0.8,1,0.05))
mtext("TMR", side=4, line=3)
legend("topleft", legend=c("PLE","TMR"), lty=c(1,2), bty="n")

################
### Figure 2 ###
################

change.LE30 <- c(as.numeric(France["PLE",])-as.numeric(France["PLE",])[1])
change.MAD <- c(as.numeric(France["MAD",])-as.numeric(France["MAD",])[1])
change.CAL <- c(as.numeric(France["CAL",])-as.numeric(France["CAL",])[1])
change.TempoLE <- c(as.numeric(France["LE*",])-as.numeric(France["LE*",])[1])

plot(my.years, change.LE30, ylim=c(0,1.5), type="l",
     xlab="Calendar Year", ylab="Change since 2009", lwd=2, main="France")
lines(my.years, change.MAD, col="gray50", lty=2, lwd=2)
lines(my.years, change.CAL, col="gray50", lty=3, lwd=2)
lines(my.years, change.TempoLE, col="gray50", lty=1, lwd=2)
legend(2010, 1.2, legend=c("PLE","MAD","CAL","LE*"), lty=c(1,2,3,1),
       col=c("black",rep("gray50",3)), bty="n", ncol=2)


#############
#Netherlands#
#############

### folder with data for the Netherlands
setwd("s:/LETHE/Geron-Tempo/Netherlands/")
HMD.Netherlands <- read.table("fltper_1x1.txt", header=TRUE, skip=2)
HMD.Netherlands$Age <- as.numeric(as.character(HMD.Netherlands$Age))
HMD.Netherlands$Age[is.na(HMD.Netherlands$Age)] <- 110


Netherlands.output <- sapply(c(2009:2015), TMR.MAD.CAL.function, HMD.frame=HMD.Netherlands)
Netherlands <- Tempo.adj.LE.function(HMD.frame=HMD.Netherlands, Output.frame=Netherlands.output)
Netherlands

################
### Figure 1 ###
################

my.years <- 2009:2015

par(mar = c(5, 4, 4, 4) + 0.1)
plot(my.years, Netherlands["PLE",], ylim=c(52.5,56.5),
     type="l", ylab="PLE", xlab="Calendar Year", main="Netherlands")
par(new = TRUE)
plot(my.years, Netherlands["TMR",], ylim=c(0.75,1), type
     = "l", axes = FALSE, bty = "n", xlab = "", ylab = "", lty=2)
axis(side=4, at = seq(0.8,1,0.05))
mtext("TMR", side=4, line=3)
legend("topleft", legend=c("PLE","TMR"), lty=c(1,2), bty="n")

################
### Figure 2 ###
################

change.LE30 <- c(as.numeric(Netherlands["PLE",])-as.numeric(Netherlands["PLE",])[1])
change.MAD <- c(as.numeric(Netherlands["MAD",])-as.numeric(Netherlands["MAD",])[1])
change.CAL <- c(as.numeric(Netherlands["CAL",])-as.numeric(Netherlands["CAL",])[1])
change.TempoLE <- c(as.numeric(Netherlands["LE*",])-as.numeric(Netherlands["LE*",])[1])

plot(my.years, change.LE30, ylim=c(0,1.5), type="l",
     xlab="Calendar Year", ylab="Change since 2009", lwd=2, main="Netherlands")
lines(my.years, change.MAD, col="gray50", lty=2, lwd=2)
lines(my.years, change.CAL, col="gray50", lty=3, lwd=2)
lines(my.years, change.TempoLE, col="gray50", lty=1, lwd=2)
legend(2010, 1.2, legend=c("PLE","MAD","CAL","LE*"), lty=c(1,2,3,1),
       col=c("black",rep("gray50",3)), bty="n", ncol=2)


#############
###The UK####
#############

### folder with data for the UK
setwd("s:/LETHE/Geron-Tempo/UK/")
HMD.UK <- read.table("fltper_1x1.txt", header=TRUE, skip=2)
HMD.UK$Age <- as.numeric(as.character(HMD.UK$Age))
HMD.UK$Age[is.na(HMD.UK$Age)] <- 110


UK.output <- sapply(c(2009:2015), TMR.MAD.CAL.function, HMD.frame=HMD.UK)
UK <- Tempo.adj.LE.function(HMD.frame=HMD.UK, Output.frame=UK.output)
UK

################
### Figure 1 ###
################

my.years <- 2009:2015

par(mar = c(5, 4, 4, 4) + 0.1)
plot(my.years, UK["PLE",], ylim=c(52.5,56.5),
     type="l", ylab="PLE", xlab="Calendar Year", main="United Kingdom")
par(new = TRUE)
plot(my.years, UK["TMR",], ylim=c(0.75,1), type
     = "l", axes = FALSE, bty = "n", xlab = "", ylab = "", lty=2)
axis(side=4, at = seq(0.8,1,0.05))
mtext("TMR", side=4, line=3)
legend("topleft", legend=c("PLE","TMR"), lty=c(1,2), bty="n")

################
### Figure 2 ###
################

change.LE30 <- c(as.numeric(UK["PLE",])-as.numeric(UK["PLE",])[1])
change.MAD <- c(as.numeric(UK["MAD",])-as.numeric(UK["MAD",])[1])
change.CAL <- c(as.numeric(UK["CAL",])-as.numeric(UK["CAL",])[1])
change.TempoLE <- c(as.numeric(UK["LE*",])-as.numeric(UK["LE*",])[1])

plot(my.years, change.LE30, ylim=c(0,1.5), type="l",
     xlab="Calendar Year", ylab="Change since 2009", lwd=2,
     main="United Kingdom")
lines(my.years, change.MAD, col="gray50", lty=2, lwd=2)
lines(my.years, change.CAL, col="gray50", lty=3, lwd=2)
lines(my.years, change.TempoLE, col="gray50", lty=1, lwd=2)
legend(2010, 1.2, legend=c("PLE","MAD","CAL","LE*"), lty=c(1,2,3,1),
       col=c("black",rep("gray50",3)), bty="n", ncol=2)

##########################
####### End of code ######
##########################
