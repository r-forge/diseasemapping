#' @rdname hook_plot_htmlsubfig
#' @export
hook_plot_array  = function (x, options) {
    # fig.tabset is a vector
    base = knitr::opts_knit$get("base.url")
    if (is.null(base)) 
        base = ""
    cap = options$fig.cap
    scap = options$fig.subcap
    
    if (is.null(cap)) 
        cap = ""
    if (is.null(scap)) {
        scap = cap
    }

    fig.num = options$fig.num
    result = sprintf("![%s](%s%s) ", scap, base, x)

    fig.names= options$fig.names
        
    fig.array = options$fig.array
    if(is.null(fig.array)) fig.array = fig.num
    Ntabset = length(fig.array)

    fig.cur = options$fig.cur
    if (is.null(fig.cur)) fig.cur = 1
    fig.index = drop(arrayInd(fig.cur, rev(fig.array)))

    theOnes = fig.index == 1
    theOnes[length(theOnes)] = FALSE
    firstNotOne = min(which(!theOnes))
    SnewHeaders = seq(from = Ntabset-1, by=-1, len = firstNotOne-1)
    toAddHere = ''
    for(D in SnewHeaders) {
      toAddHere = paste0('\n',strrep('#', D+1), ' ',  
                         fig.names[[D]][fig.index[Ntabset+1-D]], 
                         ' {.tabset}\n', toAddHere)
    }
    result = paste0('\n', toAddHere, '\n', 
                    strrep("#", Ntabset+1),  ' ', scap, ' \n', result)

    if (length(options$out.width)) {
            result = paste0(gsub("[[:space:]]*$", "", result), 
                "{width=", options$out.width, "}")
    }
    
    if(fig.cur == 1) {
      result = paste0(
        '# ', cap, ' {.tabset}', 
        result
      )
    }
    result
}