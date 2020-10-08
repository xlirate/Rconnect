// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// beyer_stretch
NumericMatrix beyer_stretch(const NumericMatrix& data, const NumericMatrix& mask, const double beta, const double z);
RcppExport SEXP _Rconnect_beyer_stretch(SEXP dataSEXP, SEXP maskSEXP, SEXP betaSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type mask(maskSEXP);
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(beyer_stretch(data, mask, beta, z));
    return rcpp_result_gen;
END_RCPP
}
// beyer_wrap
NumericMatrix beyer_wrap(const NumericMatrix& data, const NumericMatrix& mask, const double beta, const double z);
RcppExport SEXP _Rconnect_beyer_wrap(SEXP dataSEXP, SEXP maskSEXP, SEXP betaSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type mask(maskSEXP);
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(beyer_wrap(data, mask, beta, z));
    return rcpp_result_gen;
END_RCPP
}
// beyer_refect
NumericMatrix beyer_refect(const NumericMatrix& data, const NumericMatrix& mask, const double beta, const double z);
RcppExport SEXP _Rconnect_beyer_refect(SEXP dataSEXP, SEXP maskSEXP, SEXP betaSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type mask(maskSEXP);
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(beyer_refect(data, mask, beta, z));
    return rcpp_result_gen;
END_RCPP
}
// beyer_zero
NumericMatrix beyer_zero(const NumericMatrix& data, const NumericMatrix& mask, const double beta, const double z);
RcppExport SEXP _Rconnect_beyer_zero(SEXP dataSEXP, SEXP maskSEXP, SEXP betaSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type mask(maskSEXP);
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(beyer_zero(data, mask, beta, z));
    return rcpp_result_gen;
END_RCPP
}
// beyer_shrink
NumericMatrix beyer_shrink(const NumericMatrix& data, const NumericMatrix& mask, const double beta, const double z);
RcppExport SEXP _Rconnect_beyer_shrink(SEXP dataSEXP, SEXP maskSEXP, SEXP betaSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type mask(maskSEXP);
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(beyer_shrink(data, mask, beta, z));
    return rcpp_result_gen;
END_RCPP
}
// convolve_stretch
NumericMatrix convolve_stretch(NumericMatrix data, std::vector<NumericMatrix> kernel);
RcppExport SEXP _Rconnect_convolve_stretch(SEXP dataSEXP, SEXP kernelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< std::vector<NumericMatrix> >::type kernel(kernelSEXP);
    rcpp_result_gen = Rcpp::wrap(convolve_stretch(data, kernel));
    return rcpp_result_gen;
END_RCPP
}
// convolve_wrap
NumericMatrix convolve_wrap(NumericMatrix data, const std::vector<NumericMatrix>& kernel);
RcppExport SEXP _Rconnect_convolve_wrap(SEXP dataSEXP, SEXP kernelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const std::vector<NumericMatrix>& >::type kernel(kernelSEXP);
    rcpp_result_gen = Rcpp::wrap(convolve_wrap(data, kernel));
    return rcpp_result_gen;
END_RCPP
}
// convolve_refect
NumericMatrix convolve_refect(NumericMatrix data, const std::vector<NumericMatrix>& kernel);
RcppExport SEXP _Rconnect_convolve_refect(SEXP dataSEXP, SEXP kernelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const std::vector<NumericMatrix>& >::type kernel(kernelSEXP);
    rcpp_result_gen = Rcpp::wrap(convolve_refect(data, kernel));
    return rcpp_result_gen;
END_RCPP
}
// convolve_zero
NumericMatrix convolve_zero(NumericMatrix data, const std::vector<NumericMatrix>& kernel);
RcppExport SEXP _Rconnect_convolve_zero(SEXP dataSEXP, SEXP kernelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const std::vector<NumericMatrix>& >::type kernel(kernelSEXP);
    rcpp_result_gen = Rcpp::wrap(convolve_zero(data, kernel));
    return rcpp_result_gen;
END_RCPP
}
// convolve_shrink
NumericMatrix convolve_shrink(NumericMatrix data, const std::vector<NumericMatrix>& kernel);
RcppExport SEXP _Rconnect_convolve_shrink(SEXP dataSEXP, SEXP kernelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const std::vector<NumericMatrix>& >::type kernel(kernelSEXP);
    rcpp_result_gen = Rcpp::wrap(convolve_shrink(data, kernel));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Rconnect_beyer_stretch", (DL_FUNC) &_Rconnect_beyer_stretch, 4},
    {"_Rconnect_beyer_wrap", (DL_FUNC) &_Rconnect_beyer_wrap, 4},
    {"_Rconnect_beyer_refect", (DL_FUNC) &_Rconnect_beyer_refect, 4},
    {"_Rconnect_beyer_zero", (DL_FUNC) &_Rconnect_beyer_zero, 4},
    {"_Rconnect_beyer_shrink", (DL_FUNC) &_Rconnect_beyer_shrink, 4},
    {"_Rconnect_convolve_stretch", (DL_FUNC) &_Rconnect_convolve_stretch, 2},
    {"_Rconnect_convolve_wrap", (DL_FUNC) &_Rconnect_convolve_wrap, 2},
    {"_Rconnect_convolve_refect", (DL_FUNC) &_Rconnect_convolve_refect, 2},
    {"_Rconnect_convolve_zero", (DL_FUNC) &_Rconnect_convolve_zero, 2},
    {"_Rconnect_convolve_shrink", (DL_FUNC) &_Rconnect_convolve_shrink, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_Rconnect(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
