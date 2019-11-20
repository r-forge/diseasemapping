


generate_row_block_information = function(x, shared_mem_size = 1024L, showWarnings=FALSE) {

  
if(!all(class(x) %in% c('dCHMsimpl','dsRMatrix'))) {
  warning("x should be row-oriented dsRMatrix or dCHMsimpl")
  x = as(x, 'dsRMatrix')
}
  
# 'http://viennacl.sourceforge.net/doc/classviennacl_1_1compressed__matrix.html' #set
# http://viennacl.sourceforge.net/doc/compressed__matrix_8hpp_source.html, line 794

# number of column indices loaded to shared memory, number of floating point values loaded to shared memory

row_blocks_out = rep(0L, ncol(x))
row_blocks_overflow_out = rep(0L, ncol(x))
row_block_num = 1L # in R index from 1

row_blocks_out[1] = 1L
num_entries_in_current_batch = 0L;

D = 1
while(D <= nrow(x)) {
	entries_in_row = x@p[D+1] - x@p[D]
	if(entries_in_row > shared_mem_size) {
		if(showWarnings) warning('too many entries for shared memory')
		if(num_entries_in_current_batch > 0) {

			# previous rows into a block
			row_blocks_out[row_block_num] = as.integer(D-1) 
			row_blocks_overflow_out[row_block_num]=0L
			num_entries_in_current_batch = 0L

		}

		row_block_num = row_block_num + 1

		# this row in it's own block
		row_blocks_out[row_block_num] = as.integer(D) 
		# start of overflow
		row_blocks_overflow_out[row_block_num]=as.integer(shared_mem_size)
		num_entries_in_current_batch = 0L

	} else {

 	num_entries_in_current_batch = num_entries_in_current_batch + entries_in_row

 	# if this row is too big to go in current block
	if (num_entries_in_current_batch > shared_mem_size) {

		# number of rows in block before adding row D
		rows_in_batch = D - row_blocks_out[row_block_num];

		if (rows_in_batch > 0) {
			# at least one full row is in the batch. Use current row in next batch.
			D = D - 1L
			row_block_num = row_block_num + 1L
			row_blocks_out[row_block_num] = as.integer(D) 
		} else {
			# row is larger than buffer in shared memory
			row_block_num = row_block_num + 1		
			row_blocks_out[row_block_num] = as.integer(D + 1)
		}
		row_blocks_overflow_out[row_block_num]=0L
		num_entries_in_current_batch = 0L;

	} # if entries > shared mem
}
	D = D + 1

} # while loop
if (num_entries_in_current_batch > 0) {
	row_block_num = row_block_num + 1		
	row_blocks_out[row_block_num] = ncol(x)

	if (num_entries_in_current_batch < shared_mem_size) {
		row_blocks_overflow_out[row_block_num]=as.integer(0)
	} else  {
		row_blocks_overflow_out[row_block_num]=as.integer(shared_mem_size)
	}
}



list(
	row_block_num=row_block_num,
	row_blocks = row_blocks_out,
	row_blocks_overflow=row_blocks_overflow_out)

}

#' @title Sparse matrix on GPU
#'
#' @export
#' 
getVclSparseMatrix = function(x) {

  if(!all(class(x) %in% c('dCHMsimpl','dsRMatrix'))) {
      warning("x should be row-oriented dsRMatrix or dCHMSimpl")
      x = as(x, 'dsRMatrix')
    }
      
blockInfo = generate_row_block_information(x)

if(class(x) == 'dCHMsimpl') {
  col_buffer = vclVector(x@i, type = 'integer')
} else {
  col_buffer = vclVector(x@j, type = 'integer')
}

list(
	row_jumper = vclVector(x@p, type = 'integer'),
	col_buffer = col_buffer,
	elements = vclVector(x@x, type = 'double'),
	rows = x@Dim[2],
	cols = x@Dim[1],
	nonzeros = length(x@x),
	row_block_num = blockInfo$row_block_num,
	row_blocks = vclVector(as.integer(blockInfo$row_blocks), type='integer'),
	row_blocks_overflow = vclVector(as.integer(blockInfo$row_blocks_overflow), type='integer')
	)
}

#' @title Sparse multiple matrices on GPU
#'
#' @export
#' 

getVclSparseMultiMatrix = function(x, elements) {
  
  if(!all(class(x) %in% c('dCHMsimpl','dsRMatrix'))) {
    warning("x should be row-oriented dsRMatrix or dCHMSimpl")
    x = as(x, 'dsRMatrix')
  }


    if(is.vector(elements)) {
    if(length(elements)==1) {
      elements= matrix(x@x, nrow=length(x@x), ncol=elements[1])
    } else {
      diagElements = x@p[-length(x@p)]+1
      addToDiag = elements
      elements= matrix(x@x, nrow=length(x@x), ncol=length(addToDiag))
      elements[diagElements,] = elements[diagElements,] +
        matrix(addToDiag, ncol=ncol(elements), 
                 nrow = length(diagElements), byrow=TRUE)
    }
  }
  
  if(class(x) == 'dCHMsimpl') {
    col_buffer = vclVector(x@i, type = 'integer')
  } else {
    col_buffer = vclVector(x@j, type = 'integer')
  }
  
  blockInfo = generate_row_block_information(x)
  list(
    row_jumper = vclVector(x@p, type = 'integer'),
    col_buffer =  col_buffer,
    elements = vclMatrix(elements, type = gpuR::typeof(elements)),
    rows = x@Dim[2],
    cols = x@Dim[1],
    matrices = ncol(elements),
    nonzeros = nrow(elements),
    row_block_num = blockInfo$row_block_num,
    row_blocks = vclVector(as.integer(blockInfo$row_blocks), type='integer'),
    row_blocks_overflow = vclVector(as.integer(blockInfo$row_blocks_overflow), type='integer')
  )
}
