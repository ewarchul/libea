source("csa_cma_es.R")
source("ppmf_cma_es.R")



sphere_fn = function(x) { sum(x^2) }
dim = 40
x0 = rep(95, dim)
set.seed(1)

start_time <- Sys.time()
result = rb_ipop_cma_esr_ppmf(x0, sphere_fn)
print("best fiotness")
print(result$best.fitness)
end_time <- Sys.time()

print(end_time - start_time)

