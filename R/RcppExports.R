# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

.convolve_stretch <- function(data, kernel) {
    .Call(`_Rconnect_convolve_stretch`, data, kernel)
}

.convolve_wrap <- function(data, kernel) {
    .Call(`_Rconnect_convolve_wrap`, data, kernel)
}

.convolve_reflect <- function(data, kernel) {
    .Call(`_Rconnect_convolve_reflect`, data, kernel)
}

.convolve_zero <- function(data, kernel) {
    .Call(`_Rconnect_convolve_zero`, data, kernel)
}

.convolve_nan <- function(data, kernel) {
    .Call(`_Rconnect_convolve_nan`, data, kernel)
}

.convolve_shrink <- function(data, kernel) {
    .Call(`_Rconnect_convolve_shrink`, data, kernel)
}

.powered_convolve_stretch <- function(data, kernel, power = 1.0) {
    .Call(`_Rconnect_powered_convolve_stretch`, data, kernel, power)
}

.powered_convolve_wrap <- function(data, kernel, power = 1.0) {
    .Call(`_Rconnect_powered_convolve_wrap`, data, kernel, power)
}

.powered_convolve_refect <- function(data, kernel, power = 1.0) {
    .Call(`_Rconnect_powered_convolve_refect`, data, kernel, power)
}

.powered_convolve_zero <- function(data, kernel, power = 1.0) {
    .Call(`_Rconnect_powered_convolve_zero`, data, kernel, power)
}

.powered_convolve_nan <- function(data, kernel, power = 1.0) {
    .Call(`_Rconnect_powered_convolve_nan`, data, kernel, power)
}

.powered_convolve_shrink <- function(data, kernel, power = 1.0) {
    .Call(`_Rconnect_powered_convolve_shrink`, data, kernel, power)
}

.cache_samc <- function(kernel, permiability, death_rate) {
    .Call(`_Rconnect_cache_samc`, kernel, permiability, death_rate)
}

.samc_cache_sizes <- function(ca) {
    .Call(`_Rconnect_samc_cache_sizes`, ca)
}

.samc_print_cache <- function(ca) {
    invisible(.Call(`_Rconnect_samc_print_cache`, ca))
}

.samc_one_step <- function(ca, pop_in, dead_in) {
    .Call(`_Rconnect_samc_one_step`, ca, pop_in, dead_in)
}

