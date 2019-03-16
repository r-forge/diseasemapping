
#include "clRNG.hpp"
#include <clRNG/mrg31k3p.h>
#include <CL/mrg31k3pkernelString.hpp>

using namespace Rcpp;
using namespace viennacl;	
using namespace viennacl::linalg;

void convertMatclRng(Rcpp::IntegerMatrix Sin, clrngMrg31k3pStream* streams){

	int Ditem,Delement,Dcis,Dg;
	int numWorkItems = Sin.ncol();

	for(Ditem =0;Ditem < numWorkItems;Ditem++){
		for(Delement=0;Delement < 3;Delement++){

			Dcis=0;
			Dg=0;
			streams[Ditem].current.g1[Delement] = Sin(Ditem,
				Dcis*6 + Dg*3 + Delement
			);
			Dg=1;
			streams[Ditem].current.g2[Delement] = Sin(Ditem,
				Dcis*6 + Dg*3 + Delement
			);

			Dcis=1;
			Dg=0;
			streams[Ditem].initial.g1[Delement] = Sin(Ditem,
				Dcis*6 + Dg*3 + Delement
			);
			Dg=1;
			streams[Ditem].initial.g2[Delement]=Sin(Ditem,
				Dcis*6 + Dg*3 + Delement
			);

			Dcis=2;
			Dg=0;
			streams[Ditem].substream.g1[Delement]=Sin(Ditem,
				Dcis*9 + Dg*3 + Delement
			);
			Dg=1;
			streams[Ditem].substream.g2[Delement] = Sin(Ditem,
				Dcis*9 + Dg*3 + Delement
			);
		}
	}


}



void runifGpu(
	viennacl::vector_base<double> &x,
    Rcpp::IntegerMatrix streamsR,
	int numWorkItems,
	int numLocalItems,
	int ctx_id){

viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));

// create streams
size_t streamBufferSize;
clrngStatus err;
clrngMrg31k3pStream* streams = clrngMrg31k3pAllocStreams(numWorkItems, &streamBufferSize, &err);
  
convertMatclRng(streamsR, streams);

// TO DO return streams to R after calling the kernel
// transfer streams to opencl as matrix
// convert to crngMgr31k3pStream in opencl
// Create buffer to transfer streams to the device.
// change clrngMrg31k3pCopyOverStreamsFromGlobal
viennacl::vector<char> bufIn(
	ctx.create_memory_without_smart_handle(
		CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
		streamBufferSize, (void *) streams), 
	1);

// copy clrngMgr31kp to the gpu
//viennacl::vector<cl_uint> clrngMrg31k3p(18);
//copy(clrngMrg31k3p_All, clrngMrg31k3p);
//viennacl::ocl::local_mem clrngMrg31k3pLocal(clrngMrg31k3p.size()*sizeof(cl_uint));


// kernel
viennacl::ocl::program & my_prog = ctx.add_program(mrg31k3pkernelString, "my_kernel");
viennacl::ocl::kernel 
		&randomUniform = my_prog.get_kernel("mrg31k3pDoubleUint");

randomUniform.global_work_size(0, numWorkItems);
randomUniform.local_work_size(0, numLocalItems);

int Nsim = x.size();
viennacl::ocl::enqueue(randomUniform(bufIn, x, //clrngMrg31k3p, clrngMrg31k3pLocal, 
	Nsim) );

}

void runifGpuHost(
	viennacl::vector_base<double> &x){

	int D, N=x.size();
	clrngMrg31k3pStream* stream = clrngMrg31k3pCreateStreams(NULL, 1, NULL, NULL);

	for(D =0; D<N; D++) {
		x(D) = clrngMrg31k3pRandomU01(stream);
	}
}

//[[Rcpp::export]]
Rcpp::IntegerMatrix  cpp_mrg31k3pCreateStreams(
	int numWorkItems
	) {

	Rcpp::IntegerMatrix result=Rcpp::IntegerMatrix(numWorkItems,18L);

    colnames(result) = CharacterVector::create(
    	"current.g1.1", "current.g1.2", "current.g1.3", "current.g2.1", "current.g2.2", "current.g2.3",
    	"initial.g1.1", "initial.g1.2", "initial.g1.3", "initial.g2.1", "initial.g2.2", "initial.g2.3",
    	"substream.g1.1", "substream.g1.2", "substream.g1.3", "substream.g2.1", "substream.g2.2", "substream.g2.3");

	size_t streamBufferSize;
	clrngStatus err;

	int Ditem,Delement,Dcis,Dg;

	clrngMrg31k3pStream* streams = clrngMrg31k3pCreateStreams(
		NULL, numWorkItems,
		&streamBufferSize, &err);

	for(Ditem =0;Ditem < numWorkItems;Ditem++){
		for(Delement=0;Delement < 3;Delement++){

			Dcis=0;
			Dg=0;
			result(Ditem,
				Dcis*9 + Dg*3 + Delement
			) = streams[Ditem].current.g1[Delement];
			Dg=1;
			result(Ditem,
				Dcis*9 + Dg*3 + Delement
			) = streams[Ditem].current.g2[Delement];

			Dcis=1;
			Dg=0;
			result(Ditem,
				Dcis*6 + Dg*3 + Delement
			) = streams[Ditem].initial.g1[Delement];
			Dg=1;
			result(Ditem,
				Dcis*6 + Dg*3 + Delement
			) = streams[Ditem].initial.g2[Delement];

			Dcis=2;
			Dg=0;

			result(Ditem,
				Dcis*6 + Dg*3 + Delement
			) = streams[Ditem].substream.g1[Delement];
			Dg=1;
			result(Ditem,
				Dcis*6 + Dg*3 + Delement
			) = streams[Ditem].substream.g2[Delement];

		}
	}

	return result;
}

//[[Rcpp::export]]
SEXP cpp_runifGpu(
    Rcpp::S4            xR,   //vector
    Rcpp::IntegerMatrix streamsR,   //vector
    int max_global_size,     
    int max_local_size) 
{

  // data
  const bool BisVCL=1;
  const int ctx_id = INTEGER(xR.slot(".context_index"))[0]-1;


  std::shared_ptr<viennacl::vector_base<double> > x = getVCLVecptr<double>(xR.slot("address"), BisVCL, ctx_id);
  
  runifGpu(*x, streamsR, max_global_size, max_local_size, ctx_id);
 
  return(Rcpp::wrap(1.0));	
}




