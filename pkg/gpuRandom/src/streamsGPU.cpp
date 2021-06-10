

#include <clRNG/mrg31k3p.h>

#include "private.h"
#include <stdlib.h>

#if defined ( WIN32 )
#define __func__ __FUNCTION__
#endif

struct clrngMrg31k3pStreamCreator_ {
  clrngMrg31k3pStreamState initialState;
  clrngMrg31k3pStreamState nextState;
  /*! @brief Jump matrices for advancing the initial seed of streams
   */
  cl_uint nuA1[3][3];
  cl_uint nuA2[3][3];
};

#define MODULAR_NUMBER_TYPE cl_uint
#define MODULAR_FIXED_SIZE 3
#include "./modularHost.c.h"

// code that is common to host and device
#include "../include/clRNG/private/mrg31k3p.c.h"



/*! @brief Matrices to advance to the next state
 */
static cl_uint mrg31k3p_A1p0[3][3] = {
  {0, 4194304, 129},
  {1, 0, 0},
  {0, 1, 0}
};

static cl_uint mrg31k3p_A2p0[3][3] = {
  {32768, 0, 32769},
  {1, 0, 0},
  {0, 1, 0}
};


/*! @brief Inverse of mrg31k3p_A1p0 mod mrg31k3p_M1
 *
 *  Matrices to go back to the previous state.
 */
static cl_uint invA1[3][3] = {
  { 0, 1, 0 },
  { 0, 0, 1 },
  { 1531538725, 0, 915561289 }
};

// inverse of mrg31k3p_A2p0 mod mrg31k3p_M2
static cl_uint invA2[3][3] = {
  { 0, 1, 0 },
  { 0, 0, 1 },
  { 252696625, 252696624, 0 }
};


/*! @brief Default initial seed of the first stream
 */
#define BASE_CREATOR_STATE { { 12345, 12345, 12345 }, { 12345, 12345, 12345 } }
/*! @brief Jump matrices for \f$2^{134}\f$ steps forward
 */
#define BASE_CREATOR_JUMP_MATRIX_1 {  \
{1702500920, 1849582496, 1656874625}, \
{ 828554832, 1702500920, 1512419905}, \
{1143731069,  828554832,  102237247} }
#define BASE_CREATOR_JUMP_MATRIX_2 {  \
{ 796789021, 1464208080,  607337906}, \
{1241679051, 1431130166, 1464208080}, \
{1401213391, 1178684362, 1431130166} }

/*! @brief Default stream creator (defaults to \f$2^{134}\f$ steps forward)
 *
 *  Contains the default seed and the transition matrices to jump \f$\nu\f$ steps forward;
 *  adjacent streams are spaced nu steps apart.
 *  The default is \f$nu = 2^{134}\f$.
 *  The default seed is \f$(12345,12345,12345,12345,12345,12345)\f$.
 */
static clrngMrg31k3pStreamCreator defaultStreamCreator = {
  BASE_CREATOR_STATE,
  BASE_CREATOR_STATE,
  BASE_CREATOR_JUMP_MATRIX_1,
  BASE_CREATOR_JUMP_MATRIX_2
};

/*! @brief Check the validity of a seed for MRG31k3p
 */
static clrngStatus validateSeed(const clrngMrg31k3pStreamState* seed)
{
  // Check that the seeds have valid values
  for (size_t i = 0; i < 3; ++i)
    if (seed->g1[i] >= mrg31k3p_M1)
      return clrngSetErrorString(CLRNG_INVALID_SEED, "seed.g1[%u] >= mrg31k3p_M1", i);
    
    for (size_t i = 0; i < 3; ++i)
      if (seed->g2[i] >= mrg31k3p_M2)
        return clrngSetErrorString(CLRNG_INVALID_SEED, "seed.g2[%u] >= mrg31k3p_M2", i);
      
      if (seed->g1[0] == 0 && seed->g1[1] == 0 && seed->g1[2] == 0)
        return clrngSetErrorString(CLRNG_INVALID_SEED, "seed.g1 = (0,0,0)");
      
      if (seed->g2[0] == 0 && seed->g2[1] == 0 && seed->g2[2] == 0)
        return clrngSetErrorString(CLRNG_INVALID_SEED, "seed.g2 = (0,0,0)");
      
      return CLRNG_SUCCESS;
}
















clrngMrg31k3pStream* clrngMrg31k3pAllocStreams(size_t count, size_t* bufSize, clrngStatus* err)
{
  clrngStatus err_ = CLRNG_SUCCESS;
  size_t bufSize_ = count * sizeof(clrngMrg31k3pStream);
  
  // allocate streams
  clrngMrg31k3pStream* buf = (clrngMrg31k3pStream*)malloc(bufSize_);
  
  if (buf == NULL) {
    // allocation failed
    err_ = clrngSetErrorString(CLRNG_OUT_OF_RESOURCES, "%s(): could not allocate memory for streams", __func__);
    bufSize_ = 0;
  }
  
  // set buffer size if needed
  if (bufSize != NULL)
    *bufSize = bufSize_;
  
  // set error status if needed
  if (err != NULL)
    *err = err_;
  
  return buf;
}

















clrngStatus clrngMrg31k3pCreateOverStreams(clrngMrg31k3pStreamCreator* creator, size_t count, clrngMrg31k3pStream* streams)
{
  // iterate over all individual stream buffers
  for (size_t i = 0; i < count; i++) {
    
    clrngStatus err = mrg31k3pCreateStream(creator, &streams[i]);
    
    // abort on error
    if (err != CLRNG_SUCCESS)
      return err;
  }
  
  return CLRNG_SUCCESS;
}












clrngMrg31k3pStream* clrngMrg31k3pCreateStreams(clrngMrg31k3pStreamCreator* creator, size_t count, size_t* bufSize, clrngStatus* err)
{
  clrngStatus err_;
  size_t bufSize_;
  clrngMrg31k3pStream* streams = clrngMrg31k3pAllocStreams(count, &bufSize_, &err_);
  
  if (err_ == CLRNG_SUCCESS)
    err_ = clrngMrg31k3pCreateOverStreams(creator, count, streams);
  
  if (bufSize != NULL)
    *bufSize = bufSize_;
  
  if (err != NULL)
    *err = err_;
  
  return streams;
}




