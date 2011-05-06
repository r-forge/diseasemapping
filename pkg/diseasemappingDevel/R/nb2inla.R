#
#Convert nb object into inla format
#


nb2inla <-function(file, nb)
{
	
	n<-length(nb)

	if(!file.create(file))
	{
		stop("Cannot open file")
	}

	txt<-paste(n, "\n", sep="")
	cat(txt, file=file, append = TRUE)

	for(i in 1:length(nb))
	{
                temp = nb[[i]]-1;ltemp=length(temp)
                if(ltemp==1 & temp[1]==-1) {ltemp=0;temp = ""}
		txt<-paste(c(i-1, ltemp, temp), collapse=" ")
		txt<-paste(txt, "\n",sep="")
		cat(txt, file=file, append = TRUE)
                
	}
}

