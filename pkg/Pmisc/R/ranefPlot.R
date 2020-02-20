#' Dotplot of random effects from glmmTMB
#' 
#' @description Creates a nice dotplot with CI's
#'
#' @param x model output from glmmTMB
#' @param component model component
#' @param grpvar name of random effect
#' @param cutSe only plot levels where conditional SE is above this threshold
#' @param maxNames maximum number of labels on the plot
#' @param col vector of colours or name of a colorBrewer pallet
#' @param xlab x label
#' @return nothing
#' @export


ranefPlot = function(x, 
  component=c('cond','zi'),
grpvar = 1,
term='(Intercept)',
cutSe = Inf,
maxNames = 40,
col = 'Dark2',
xlab = 'x', level = 0.95, ...) {

component = component[1]
if(is.numeric(grpvar))
	grpvar = names(glmmTMB::VarCorr(x)[[component]])[grpvar]
x= as.data.frame(glmmTMB::ranef(x))

x = x[x$component==component & x$grpvar == grpvar &
  x$term == term, ]

if(length(level) == 1) {
  level = (1-level)/2
  level = c(level, 1-level)
}


x$lower = x$condval + qnorm(level[1])*x$condsd
x$upper = x$condval + qnorm(level[2])*x$condsd
x = x[x$condsd < cutSe, ]
x = x[order(x$condval), ]

x$index = rank(x$condval)
x$accurate = rank(x$condsd) < maxNames

if(all(col %in% rownames(RColorBrewer::brewer.pal.info)))
 col = RColorBrewer::brewer.pal(8, 'Dark2')

x$col= rep_len(col, nrow(x))
x$colTrans = paste0(x$col, '40')
x$colLine = x$col
x[!x$accurate,'colLine'] = x[!x$accurate,'colTrans']

x$cex = -sqrt(x$condsd) 
x$cex = x$cex - min(x$cex) + 0.1/3
x$cex = 3*x$cex / max(x$cex)


x$textpos = rep_len(c(4,2), nrow(x))
x[!x$accurate & x$condval > 0, 'textpos'] = 4
x[!x$accurate & x$condval < 0, 'textpos'] = 2

x$textloc = x$condval


x$textCex = c(0.5, 0.9)[1+x$accurate]


par(mar=c(4,0,0,0), bty='n')

forPlot = list(x=x$condval, y=x$index, 
  yaxt='n', xlim = c(min(x[x$accurate, 'lower']), max(x[x$accurate, 'upper'])),
  xlab=xlab, ylab='', pch=15, col=x$colTrans , cex=x$cex)
theDots = list(...)
forPlot = forPlot[setdiff(names(forPlot), names(theDots))]

do.call(plot, c(forPlot, theDots))

#plot(x$condval, x$index, yaxt='n', xlim = range(x$condval)*1.2,
 #  xlab=xlab, ylab='', pch=15, col=x$colTrans , cex=x$cex, ...)

x[!x$accurate & x$condval > 0, 'textloc'] = par('usr')[1]
x[!x$accurate & x$condval < 0, 'textloc'] = par('usr')[2]

abline(v=0, col='grey')
segments(x$lower, x$index, x$upper, x$index, pch=15, col=x$colLine)
text(
  x$textloc, 
  x$index, x$grp,
  pos = x$textpos,
  col=x$col,
  cex=x$textCex, offset=1)
}
