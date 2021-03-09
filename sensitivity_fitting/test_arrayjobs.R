# lets again using freeze function to test how does the array jobs work

# first a simple array job where each subtask using only one core 

# args <- as.numeric( commandArgs( TRUE ) )
# 
# job_id <- args[1]
# print(job_id)
# library(tictoc, lib.loc = "/projappl/project_2003466/project_rpackages")
# 
# 
# tic()
# Sys.sleep(time = job_id/job_id * 10)
# toc()

# second try an array job where each sub task use 10 cores to parallel future_map frezzing

args <- as.numeric( commandArgs( TRUE ) )

job_id <- args[1]
print(job_id)
library(tictoc, lib.loc = "/projappl/project_2003466/project_rpackages")
library(furrr)

tic()
plan(multicore,workers=6)
x= c(10,10,10,10,10,10)* (job_id/job_id)
nothingness <- future_map(.x = x, .f = ~Sys.sleep(.x))
plan(sequential)
toc()
