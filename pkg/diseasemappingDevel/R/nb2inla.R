#
#Convert nb object into inla format
#Written By Virgilio Gomez-Rubio


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
		txt<-paste(c(i-1, length(nb[[i]]), nb[[i]]-1), collapse=" ")
		txt<-paste(txt, "\n",sep="")
		cat(txt, file=file, append = TRUE)
	}
}

