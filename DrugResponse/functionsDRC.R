
cellCountCV <- function(rawCVdata, cellnumcol='Cell.Nucleus', normPos=1) {
  #' For 96 well and 384 well plates
  #' Function to process cell count data in Cellavista multi-measurement format
  #'
  #' @param rawCVdata \emph{data.frame}
  #' @param cellnumcol \emph{character}
  #' @param normPos \emph{integer}
  #'
  #' Function takes \code{data.frame} of raw cell count data exported from Synentec
  #' Cellavista and extracts time of image acquisition, identifies wells from which images
  #' have been obtained, outputs a \code{data.frame} with \code{colnames} of
  #' \emph{Row, Column, Well, Cell.count, nl2}. \emph{nl2} is log2 values of \emph{Cell.count}
  #' normalized to \emph{normPos} argument (default = 1).
  #'
  #' Default column name containing cell counts is \emph{Cell.Nucleus} which is from
  #'  the \code{Cell Nuclei v1.0} image processing function on the Cellavista instrument.
  out <- data.frame()
  
  d <- rawCVdata
  d <- d[,c('Column','Row',cellnumcol)]
  # rename Cell.Nucleus to Cell.count
  colnames(d)[3] <- 'Count'
  
  # ensure that Row and Column are char vectors
  d$Row <- as.character(d$Row)
  d$Column <- as.character(d$Column)
  
  # make row index of data positions for each time point
  idx <- grep('Timespan',d[,1])
  
  # extract time of acquisition for each time point
  times <- as.numeric(d[,2][idx])
  time.idx <- seq(length(times))
  
  # extract the number of wells we have data for
  numWells <- idx[2]-idx[1]-2
  
  # assemble data into single structure (out)
  for(tp in time.idx)
  {
    temp <- d[(idx[tp]+1):(idx[tp]+numWells),]
    temp$Time  <- times[tp]
    out <- rbind(out,temp)
  }
  
  rownames(out) <- seq(nrow(out))
  
  # Generate a column for "Well" in the format RCC
  out  <- fixCVWellName(out)
  out$Well <- paste0(out$Row,"-",out$Column)
  
  return(out)
}

fixCVWellName <- function(countData){
  #For 96 well and 384 well plates
  rowString = c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P")
  
  for(i in 1:length(rowString)){
    if(i < 10){
      countData[countData$Row == rowString[i],"Row"] = paste0("R0",i)
    } else{
      countData[countData$Row == rowString[i],"Row"] = paste0("R",i)
    }
  }
  
  for(j in 1:24){
    if(j < 10){
      countData[countData$Column == j,"Column"] = paste0("C0",j)
    } else{
      countData[countData$Column == j,"Column"] = paste0("C",j)
      }
  }
  return(countData)
}

remDrugChange <- function(countdata, interval = 4){
  #largest interval between drug changes, but smaller than interval between images
  
  timelist = unique(countdata$Time)
  discTime = list()
  j=1
  #Finds imaging times that have small intervals between them. Writes first time to a list.
  for(i in 1:(length(timelist)-1)){
    diff = timelist[i+1]-timelist[i]
    if(diff <= interval){
      discTime[j] = i
      j = j + 1
    } else{}
  }
  dTime = timelist[as.numeric(discTime)]
  #Remove the first time from the dataframe. This is the drug change time point.
  subCountDat = subset(countdata,!countdata$Time%in%dTime)
  return(subCountDat)
}


expSmooth <- function(countdata, alpha = 0.1, countcol = "Count"){
  #Use this on labeled count data before nl2 and stats
  #Alpha is the smoothing coeficcient. Lower alpha, higher dampening (1-alpha).
  #Specify countcol.
  Wells = unique(countdata$Well)
  newdf = data.frame()
  for(i in 1:length(Wells)){
    sdf = subset(countdata,countdata$Well==Wells[i])
    #order sdf by time increasing
    sdf = sdf[order(sdf$Time,decreasing = F),]
    times = unique(sdf$Time)
    #Empty list for smoothed values
    svals = c()
    for(j in 1:length(times)){
      if(j == 1){
        svals[j] = sdf[sdf$Time==times[j],countcol]
      } else{
        #Get last smoothed value
        lsval = svals[j-1]
        #Get current value
        cval = sdf[sdf$Time==times[j],countcol]
        #Make exponentially smoothed value
        svals[j] = alpha*cval + (1-alpha)*lsval
      }
    }
    #Put list of svals into dataframe and rbind for newdf
    sdf$esCount = svals
    newdf = rbind(newdf,sdf)
    rm(sdf,svals)
  }
  
  return(newdf)
  
}



countImagJ <- function(platefile){ 
  
  plate1 = read.delim(platefile, header = T, sep="\t")
  plate1$Datetime = substr(plate1$Slice, 1, 14)
  plate1$Datetime = ymd_hms(plate1$Datetime)
  times = difftime(plate1$Datetime, plate1$Datetime[1], units = "hours")
  plate1$Time = as.numeric(times)
  plate1$Well = substr(plate1$Slice, (nchar(as.character(plate1$Slice)) - 6), nchar(as.character(plate1$Slice)))
  
  cellcount = aggregate(Count~ Time + Well, data = plate1, FUN = 'sum')
  cellcount$Row    = substr(cellcount$Well, 1, 3)
  cellcount$Column = substr(cellcount$Well, 5, 7)
  
  return(cellcount)
  
}


platDrug <- function(platemap, countdata){
  
  rows = as.character(platemap$X)
  cols = colnames(platemap)[2:12]
  rownames(platemap) = rows
  platemap$X = NULL
  
  for(i in 2:length(rows)){
    for(j in 2:length(cols)){
      Well = paste0(rows[i],"-",cols[j])
      countdata[countdata$Well == Well,"Drug"] = platemap[i,"Drug"]
      countdata[countdata$Well == Well,"Conc"] = platemap[i,j]
      countdata[countdata$Well == Well,"Units"] = platemap[i,"Conc"]
    }
  }
  
  return(countdata)
}

styleDrugPlate96 <- function(platemapfile, countdata){

  #Get plate rows and columns. Make vectors.
  plateDat = read_excel(platemapfile, sheet=1, col_names = F)
  rows = pull(plateDat[2:9,1])
  # rows = as.character(plateDat[2:9,1])
  # rows = rows[!is.na(rows)]
  cols = as.character(plateDat[1,2:13])
  
  #List colors that match conditions and iterate
  formats = xlsx_formats(platemapfile)
  x = xlsx_cells(platemapfile)
  colList = formats$local$fill$patternFill$fgColor$rgb
  colList = colList[!is.na(colList)]
  
  #Label or fill plateDat$Color with color style
  for(k in 1:length(colList)){
    style = x[x$local_format_id %in% which(formats$local$fill$patternFill$fgColor$rgb == colList[k]),]
    #Take col > 12 and find plateDat row. This is where the condition matrix is on the excel file.
    style = style[style$col>12,]
    uRow = unique(style$row)
    uStyle = unique(style$local_format_id)
    #Label the color column on plateDat with the color ID
    plateDat[uRow[1],"Color"] = uStyle[1] 
  }
  
  #Assemble the condition Matrix that exists past column 13
  condMatrix = plateDat[,-c(1:13)]
  #Get intended column names
  colConMat = as.character(condMatrix[1,c(1:ncol(condMatrix)-1)])
  condMatrix = condMatrix[-1,]
  #Add column names for conditions without replacing "Color"
  colnames(condMatrix)[1:ncol(condMatrix)-1] = colConMat
  
  #Step through plate by Rxx-Cxx
  for(i in 2:length(rows)){
    for(j in 2:length(cols)){
      #Label countdata concentration listed on plateDat
      Well = paste0(rows[i],"-",cols[j])
      countdata[countdata$Well == Well,"Conc"] = plateDat[i+1,j+1]
      #Get color label of the plateDat well
      g = as.numeric(x[which(x$row==i+1 & x$col == j+1),"local_format_id"])
      #Get the conditions for the local_format_id
      condSamp = condMatrix[condMatrix$Color == g,]
      colcondSamp = colnames(condSamp)
      #Use loop to fill conditions from condMatrix to countdata
      for(f in 1:ncol(condSamp)-1){
        countdata[countdata$Well == Well, colcondSamp[f]] = condSamp[1,f]
      }
    }
  }
  
  return(countdata)
  
}


getPlateGroupString96 <- function(platemapfile){
  plateDat = read_excel(platemapfile, sheet=1, col_names = F)
  #Assemble the condition Matrix that exists past column 13
  condMatrix = plateDat[,-c(1:13)]
  #get groups/conditions from matrix
  condHead = as.character(condMatrix[1,])
  condGroup = subset(condHead, !condHead%in%"Units")
  #Add in time for the DRC conditions/groupings
  condGroup = c("Time",condGroup, "Conc")
  return(condGroup)
}


#Normalized Log2 count

compNL2Stats <- function(countdata, plateGroups, CountVar = "Count", time){
  dn = data.frame()
  cellcount = data.frame(countdata)
  cellcount$l2 = log2(cellcount[,CountVar])
  
  for(j in unique(cellcount$Well)){
    temp = subset(cellcount, cellcount$Well==j)
    temp$nl2 = temp$l2 - temp$l2[temp$Time==time]
    dn = rbind(dn, temp)
  }
  dn$Time = round(dn$Time, 4)
 
  #Compute summary stats
  drte =  summarySE(dn, measurevar="nl2", groupvars= plateGroups,
                    na.rm = T, conf.interval = 0.95, .drop = T)
  drte = drte[rowSums((is.na(drte)))==0,]
  #Plot normalized log2 count curve
  drte$Time = as.numeric(drte$Time)
  return(drte)
}



compNL2 <- function(countdata, ntimepoint = 1){
  dn = data.frame()
  cellcount = data.frame(countdata)
  cellcount$l2 = log2(cellcount$Count)
  timelist = unique(cellcount$Time)
  normTime = timelist[ntimepoint]
  for(j in unique(cellcount$Well)){
    temp = subset(cellcount, cellcount$Well==j)
    temp$nl2 = temp$l2 - temp$l2[temp$Time==normTime]
    dn = rbind(dn, temp)
  }
  dn$Time = round(dn$Time, 4)
  return(dn)
}