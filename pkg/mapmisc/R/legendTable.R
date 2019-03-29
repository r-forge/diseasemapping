legendTable = function(x,
    type=c('latex', 'html'),
    box = c(-0.2, 1, 2),
    unit = 'em',
    collapse=NULL) {
  
  type = type[1]
  
  if(length(grep("^Raster",class(x))))
    x = levels(x)[[1]]
  
  if(length(unit)>1 & type %in% names(unit)){
    unit = unit[type]
  }

  if(!length(x$label)) {
		x$label = x$legend
	}
		
  if(!length(x$label))
    x$label = paste(
      '[',
      x$breaks[-length(x$breaks)],
        ', ',
        x$breaks[-1],
    ']'
        )
  
  if(type=='latex'){
    res = legendTableLatex(x, box, unit, collapse)
  }
  if(type=='html'){
    res = legendTableHtml(x, box, unit, collapse)
  }
  res
}




legendTableHtml = function(x, box, unit, collapse) {
  
  box = box[length(box)]
  
	if(is.numeric(box)) {
	Scol = paste(
      '<span style="background-color: ',
      substr(x$col,1,7), '; ',
      'padding-right: ', box, 
      'em"></span>', sep='')
	} else {
		Scol = paste(
      	'<span style="color: ',
      	substr(x$col,1,7),  
      	'">', box, 
				'</span>', sep='')
	}

  thetable=data.frame( col = Scol)

  thetable$label = x$label
	thetable = thetable[!is.na(x$col),]

		
  if(length(collapse))
    thetable = paste(
        paste(
        		thetable$label, ' (',
        		thetable$col, ')'
        ),
      	collapse=collapse
    )
	
  thetable
}

legendTableLatex = function(x, rule, ruleUnit, collapse){
  
  if(length(rule)==1) rule = c(0, rule)
  if(length(rule)==2) rule = c(rule, rule[2]*2)
  if(is.numeric(rule)) rule = paste(rule, ruleUnit, sep='')
  
  latexCol = paste(
      '\\protect{\\textcolor[HTML]{',
      substr(x$col,2,7),
      '}{\\rule[',
      rule[1],']{',
      rule[3],'}{',
      rule[2],'}}}',
      sep=''
  )
  
  thetable=data.frame(col=latexCol, stringsAsFactors=FALSE)

  
  thetable$label = gsub("\\%", "pct", x$label)
  thetable$label = gsub("\\{|\\}", "", thetable$label)
  thetable = thetable[!is.na(x$col),]
  
  if(length(collapse))
    thetable = paste(
        paste(
        thetable$label, ' (',
        thetable$col, ')'
        ),
      collapse=collapse
    )

  thetable
}