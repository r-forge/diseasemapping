#include "geostatsgpu.hpp"


template<typename T>
double luT(
	viennacl::matrix<T> &vclX,
	viennacl::vector_base<T> &vclD) {

	double logdet=-99.9;

	viennacl::linalg::lu_factorize(vclX);

	// pointer to the actual diagonal
	viennacl::vector_base<T> diagOfVar(
			vclX.handle(), vclX.size1(), 0, 
			vclX.internal_size2() + 1);

		// compute log determinant
		vclD = viennacl::linalg::element_log(diagOfVar);
		logdet = viennacl::linalg::sum(vclD);
// OPERATION_UNARY_LOG_TYPE 	
		//http://viennacl.sourceforge.net/doc/scheduler_8cpp-example.html#a11

		// put the diagonals in D, and 1's on the diagonal of L
		vclD = diagOfVar;
		diagOfVar = 1.0;

		return(logdet);
}

template<typename T>
double cpp_luT(
	Rcpp::S4 xR,
	Rcpp::S4 dR) {

	const bool BisVCL=1;
	const int ctx_id = INTEGER(xR.slot(".context_index"))[0]-1;

	std::shared_ptr<viennacl::matrix<T> > vclX = getVCLptr<T>(xR.slot("address"), BisVCL, ctx_id);
	std::shared_ptr<viennacl::vector_base<T> > vclD = getVCLVecptr<T>(dR.slot("address"), BisVCL, ctx_id);


	return(luT(*vclX, *vclD));
}

//[[Rcpp::export]]
SEXP cpp_lu(
	Rcpp::S4 xR,
	Rcpp::S4 dR) {

	const std::string theClass =  Rcpp::as<std::string>(xR.attr("class"));
	double logDet;

    if(theClass == "fvclMatrix") {
        logDet = cpp_luT<float>(xR, dR);
    } else if(theClass == "dvclMatrix") {
        logDet = cpp_luT<double>(xR, dR);
    } else {
        Rcpp::exception("unknown type detected for x");
        logDet = -99.9;
    }
	return(Rcpp::wrap(logDet));
}