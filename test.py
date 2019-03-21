import pipeline_masse_dans_coeur as pmdc

simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire10_hr'
owner = 'averliat'
num_output = 55

rho_seuil = 1e11 #particules par cm^3

res=pmdc.masse_dans_coeur_output_unique(simu, owner, num_output, rho_seuil)

