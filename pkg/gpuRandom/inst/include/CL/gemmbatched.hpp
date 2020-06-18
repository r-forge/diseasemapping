//A is (M*z)*K
//B is (K*z)*N
//group size is M*N


// Tiled and coalesced version
__kernel void myGEMM2(const int M, 
                      const int N, 
                      const int K,
                      const int z,
                      const __global float* A,
                      const __global float* B,
                      __global float* C) {
  
  // Thread identifiers
  const int localrow = get_local_id(0); // Local row ID 
  const int localcol = get_local_id(1); // Local col ID 
  
  // Identify the row and column of the C matrix to work on
  const int globalRow = M*get_group_id(0) + localrow; // Row ID of C 
  const int globalCol = localcol; // Col ID of C 
  
  // Local memory to fit a tile of TS*TS elements of A and B
  __local float Asub[M][K];
  __local float Bsub[K][N];
  
  // Initialise the accumulation register
  float acc = 0.0f;
  
  // Loop over all batches
  //const int numBatches = z;
  for (int t=0; t<z; t++) {
    
    // Load one batch of A and B into local memory
    const int batchedRow = M*t + localrow;
    const int batchedCol = localcol
      Asub[col][row] = A[batchedCol*M*z + globalRow];
    Bsub[col][row] = B[globalCol*K*z + batchedRow];
    
    // Synchronise to make sure the tile is loaded
    barrier(CLK_LOCAL_MEM_FENCE);
    
    // Perform the computation for a single tile
    for (int k=0; k<K; k++) {
      acc += Asub[k][row] * Bsub[col][k];
    }
    
    // Synchronise before loading the next tile
    barrier(CLK_LOCAL_MEM_FENCE);
  }
  
  // Store the final result in C
  C[globalCol*M*z + globalRow] = acc;
}

























