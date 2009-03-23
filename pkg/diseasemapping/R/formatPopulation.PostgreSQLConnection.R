formatPopulation.PostgreSQLConnection <- function(dbname,user="postgres",driver="PostgreSQL",aggregate=FALSE) {

  # function to get data from a database
  # Input:
  # dbname = name of database
  # user = user name for the database, default='postgres'
  # driver = database driver, default='PostgreSQL'
  # aggregate = if TRUE, return aggregated population data by age and sex group
  # Output:
  # the raw data from the database, or the aggregated data if requested

  # functionality very limited, don't think this can work on data except our current data

  # example
  #result <- formatPopulation.PostgreSQLConnection("spatial")
  #result.agg <- formatPopulation.PostgreSQLConnection("spatial",aggregate=T)

  
  library(RPostgreSQL)

  # should get this from the database!
  stuff = c("CSDUID","CSDNAME","CSDTYPE","PRUID","PRNAME","CDUID","CDNAME","CDTYPE","CMAUID","CMANAME","SACTYPE","ERUID","ERNAME","TotPop2001","TotPop2006","TotPopPerc","TotPop","AreaSqKm20","TotPop_x_A","MaleTot","FemaleTot")
  
  sql <- c()

  begin <- TRUE # start loop indicator

  for (sex in c("M","F")) { # over Males and Females
    
    for (age in seq(0,80,by=5)) { # need to handle M85plus and F85plus!
    
      if (aggregate==TRUE) {

        # sum things up for each category
        one <- paste('select sum("',sex,age,'_',age+4,'") as population,\'',sex,'\' as sex, \'',age,'_',age+4,'\' as age from "Ontario"."CSD2006"' ,sep="")

        if (begin==TRUE) { # turn loop start indicator off
          sql <- one
          begin <- FALSE
        }
        else sql <- paste(sql,'union',one) # combine sql commands with 'union'
        
      }

      else {
    # get the columns
        one <-   paste('select',paste(c('"',paste(stuff,collapse='","'),'",'),collapse=''))

    # add sex and age columns
        two <-   paste(c('"',sex,age,'_',age+4,'" as population, \'',sex,age,'_',age+4,'\' as group, \'',sex,'\' as sex, \'',age,'_',age+4,'\' as age from "Ontario"."CSD2006"'),collapse="")

        if (begin==T) { # turn loop start indicator off
          sql <- paste(one,two)
          begin <- F
        }
        else sql <- paste(sql,'union',one,two) # combine sql commands with 'union'
        
      }
    }
  }

  # send sql query to database
  drv <- dbDriver(driver)
  con <- dbConnect(drv,dbname=dbname,user=user)
  fs <- dbSendQuery(con,sql)
  fetch(fs,n=-1)

}



#----------- for reference ------------#

#select 'M0_4', "CSDUID", "CSDNAME", "CSDTYPE", "PRUID", "M0_4" as population, 'M' as sex, 0 as age from "Ontario"."CSD2006"
#union 
#select 'M5_9', "CSDUID", "CSDNAME", "CSDTYPE", "PRUID", "M5_9" as population, 'M' as sex, 5 as age from "Ontario"."CSD2006"
#union
#select 'M10_14', "CSDUID", "CSDNAME", "CSDTYPE", "PRUID", "M10_14" as population, 'M' as sex, 10 as age from "Ontario"."CSD2006";
