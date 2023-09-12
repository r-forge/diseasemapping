
legendBreaks = function(pos,
    breaks,
    col,
    legend,
    rev=TRUE,
    outer=TRUE,
    pch=15,
    bg='white',
    cex=par('cex'),
    pt.cex=2.5*cex,
    text.col=par('fg'),
    title=NULL,
    inset=0.05,
    title.col=text.col,
    adj=0,
    width=Inf, 
    lines=Inf,
    y.intersp,
    ...){
  
  if(!missing(breaks)){
    if(is.factor(breaks)){
      # if it's a raster
      if(length(grep("^Raster",class(breaks)))){
        breaks = levels(breaks)[[1]]
      } else {
        breaks=list(legend=levels(breaks))
      }
    }
  }
  
  if( missing(legend) & missing(breaks))
    warning("legend or breaks must be supplied")
  if(missing(legend)&!missing(breaks)) {
    if(is.list(breaks)){
      legendCol = intersect(
          c('legend','label','level','breaks','ID'),
          names(breaks)
      )
      if(!length(legendCol)){
        warning("can't find legend in breaks")
      }
      legend = breaks[[ legendCol[1] ]]
    } else { # breaks isn't a list (or df)
      legend=breaks
    }
  }
  
  if(missing(col)){
    col='black'
    if(!missing(breaks)) {
      if(is.list(breaks)) {
        if(any(names(breaks)=='col'))
          col = breaks[['col']]
      }
    }
  }
  
  if(rev){
    col=rev(col)
    legend=rev(legend)
  }
  
  diffYmult = 0
  if(length(col) == (length(legend)-1)) {
    # one more legend item than colours
    col = c(NA, col)
    pch = c(NA,
        pch[round(seq(1, length(pch), len=length(legend)-1))]
    )
    diffyMult=1
    # make text transparent, add legend text manually afterwards
    # because graphics::legend doens't align it correctly
    theTextCol = '#FFFFFF00'
  } else { # same number of colours as legend entries
    theTextCol = text.col
    # get rid of entries where col is NA
    theNA = is.na(col)
    if(any(theNA)){
      col = col[!theNA]
      legend = legend[!theNA]
    }
  }
  
# line wrapping for legend labels
  if(any(nchar(as.character(legend)) > width)) {
    legend =  trimws(
        gsub(
            paste('(.{1,', width, '})(\\s|/|$)' ,sep=''), 
            '\\1\n ', 
            as.character(legend)
        )
    )
  }
  
  # remove excess lines
  theNewLines = gregexpr('\n', as.character(legend))
  toCrop = which(unlist(lapply(theNewLines, length)) >= lines)
  if(length(toCrop)) {
    cropPos = unlist(lapply(theNewLines[toCrop], function(qq) qq[lines]))
    legend = as.character(legend)
    legend[toCrop] = 
        trimws(substr(legend[toCrop], 1, cropPos))
  }
  
  shiftLegendText = rep(0, length(legend))
  
  if(missing(y.intersp)){
    
    if(is.character(legend)) {	
      theNewLines = gregexpr('\n', as.character(legend))
      y.intersp=max(
          c(1.25, 
              0.5+unlist(
                  lapply(theNewLines, function(qq) sum(qq>0))
              )
          ) 
      ) - 0.25
    } else {
      y.intersp = 1
      
      # format legend as character
      # note, if y.intersp is supplied this won't be done
      if(is.numeric(legend)) {
        
        legend = as.character(legend)
        widthHere = strwidth(legend, cex=cex)
        maxWidth = max(widthHere) 

        # padd for minus sign
        withMinus = grep("^[[:space:]]*[-]", legend)
        toAddForMinus = rep(0, length(legend))
        toAddForMinus[-withMinus] = pmin(
            maxWidth - widthHere[-withMinus],
            strwidth("-", cex=cex)
        )
        # width before decimal
        charNoDec = strwidth(gsub("(e|[.])[[:digit:]]*$", "", legend), cex=cex)
        maxCharNoDec = max(charNoDec)
        toAddLeft = pmin(
            maxCharNoDec - charNoDec,
            maxWidth - widthHere)
        
        # width after decimal
        Ndec = strwidth(
            gsub("^[[:space:]]*[[:punct:]]*([[:digit:]]|e[+])+ *", "", legend), 
            cex=cex)
        maxDec = max(Ndec)
        
        toAddRight = pmin(
            maxDec - Ndec,
            maxWidth - widthHere)

        # ideally we'd add space for minus, and padding before decimal, and padding after decimal
        idealWidth = widthHere  + toAddRight + toAddLeft #+ toAddForMinus
        tooWide = idealWidth - maxWidth
        
        shiftLegendText  = pmin(
            toAddForMinus + toAddLeft - 0.4*tooWide,
            maxWidth - widthHere   
        )
            


        
      } # end justification
      
    }
  }
  if(all(is.na(y.intersp))){
    y.intersp=1
  }
  adj = rep_len(adj, 2)
  adj[2] = adj[2] + y.intersp/4
  
  # get rid of transparency in col
  withTrans = grep("^#[[:xdigit:]]{8}$", col)
  col[withTrans] = gsub("[[:xdigit:]]{2}$", "", col[withTrans])
  
  if(outer){
    oldxpd = par("xpd")
    oldpar = par(no.readonly = TRUE) 
    par(xpd=NA)
  on.exit(par(oldpar))        
  fromEdge = matrix(par("plt"), 2, 2, 
        dimnames=list(c("min","max"), c("x","y")))
    propIn = apply(fromEdge, 2, diff)
    if(is.character(pos)) {
      forInset = c(0,0)
      if(length(grep("left$", pos))){
        forInset[1] = -fromEdge["min","x"]					
      } else if(length(grep("right$", pos))){
        forInset[1] = fromEdge["max","x"]-1					
      }
      if(length(grep("^top", pos))){
        forInset[2] = -fromEdge["min","y"]					
      } else if(length(grep("^bottom", pos))){
        forInset[2] = fromEdge["max","y"]-1					
      }
      
      inset = forInset/propIn + inset
    }
  }
  
#	legend = format(as.character(legend), justify='right')
  result=legend(
      pos,
      legend=legend,
      bg=bg,
      col=col,
      pch=pch,
      pt.cex=pt.cex,
      inset=inset,
      cex=cex,
      text.col=theTextCol,
      title.col=title.col,
      title=title,
      y.intersp=y.intersp,
      adj=adj,
      ...
  )
  
  if(text.col != theTextCol) {
    diffy = diff(result$text$y)/2
    diffy = c(
        diffy,diffy[length(diffy)]
    )*diffyMult
    result$text$y = result$text$y + diffy
    
    
    result$text$xOrig = result$text$x 
    result$text$x = result$text$x + shiftLegendText/2 + max(strwidth(legend, cex=cex))/2
    
    if(par("xlog")) result$text$x = 10^result$text$x
    if(par("ylog")) result$text$y = 10^result$text$y
    
    
    graphics::text(
        result$text$x, 
        result$text$y,
        legend, 
        col=text.col,
        adj=0.5,
        cex=cex)
  }      
  
  result$legend = legend
  return(invisible(result))
}
