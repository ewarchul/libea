source("csa_cma_es.R")



sphere_fn = function(x) { sum(x^2) }
dim = 1000
x0 = rep(95, dim)

start_time <- Sys.time()
result = rb_ipop_cma_esr_csa(x0, sphere_fn)
end_time <- Sys.time()

print(end_time - start_time)

