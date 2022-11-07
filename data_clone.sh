
# first load allas module and configure according to the project folder name:
#https://docs.csc.fi/data/Allas/accessing_allas/

# module load allas
# allas-conf project_2003466

cd /scratch/project_2003466
rclone copy ./ces_21q1_io/ allas:drug_moa/ces_21q1_io/
rclone copy ./depmap_21q1/ allas:drug_moa/depmap_21q1/
rclone copy ./depmap_static/ allas:drug_moa/depmap_static/

rclone copy ./ctrp/ allas:drug_moa/ctrp/
rclone copy ./gdsc/ allas:drug_moa/gdsc/
rclone copy ./prism/ allas:drug_moa/prism/
rclone copy ./forward_modelling/ allas:drug_moa/forward_modelling/
rclone copy ./L1000/consensus_signatures_challenge/ allas:drug_moa/L1000/consensus_signatures_challenge/
rclone copy ./prior/ allas:drug_moa/prior/
rclone copy ./unsupervised_target_pred/ allas:drug_moa/unsupervised_target_pred/
rclone copy ./glmnet_modelling/target_pred_server/res_pathway_all.RData allas:drug_moa/unsupervised_target_pred/res_pathway_all.RData

