library(tidyverse)
library(tidymodels)

load("~/cluster_scratch/forward_modelling/gdsc_modelres/drug_1.RData")
# save.image("~/cluster_scratch/forward_modelling/forwardmodelling_all_update.RData")
ces1_perf_new <- ces1_perf 
load("~/cluster_scratch/forward_modelling/gdsc_spearman/drug_1.RData")
ces1_perf_old <- ces1_perf 


ces1_perf_new <- ces1_perf_new %>% 
  mutate(spearman_cor= map_dbl(.x = final_res,
                               .f =function(x){cor(x$y, x$.pred,method= "spearman")}
      )
  )
    
    

load("~/cluster_wrk/drug_moa/sensitivity_fitting/forward_modelling/tmp.RData")
# different result is because augment should not be used to glmnet object
# it allows glm object but not elnet object?
tmp4 <- tmp2 %>% pull_workflow_fit()
predict(tmp2, new_data = tmp3)
tmp5 <- bake(object = tmp2 %>% pull_workflow_prepped_recipe(),new_data= tmp3)

predict(tmp4, new_data= tmp5)
predict(tmp2, new_data= tmp3)



# https://broom.tidymodels.org/reference/tidy.glmnet.html

# https://broom.tidymodels.org/reference/augment.lm.html