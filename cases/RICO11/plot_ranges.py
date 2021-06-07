xscaledict = {
  "thl" : "linear",
  "00rtot" : "linear",
  "rliq" : "linear",
  "prflux" : "linear",
  "cl_nc" : "linear",
  "cl_nr" : "linear",
  "clfrac" : "linear",
  "wvar" : "linear",
  "w3rd" : "linear",
  "sat_RH" : "linear",
  "rad_flx" : "linear",
  "lwp" : "linear",
  "rwp" : "linear",
  "er" : "linear",
  "wvarmax" : "linear",
  "surf_precip" : "linear",
  "acc_precip" : "linear",
  "cloud_base" : "linear",
  "gccn_rw_cl" : "linear",
  "non_gccn_rw_cl" : "linear",
  "base_prflux_vs_clhght" : "log",
  "base_prflux_vs_clhght number of occurances" : "log",
  "cl_gccn_conc" : "linear",
  "clb_bigrain_mean_rd" : "linear", 
  "clb_bigrain_mean_kappa" :  "linear", 
  "clb_bigrain_mean_conc" :  "linear", 
  "clb_bigrain_mean_inclt" :  "linear", 
  "clb_bigrain_mean_gccn_fraction" :  "linear", 
  "cloud_cover_rico" : "linear",
  "min_cloud_base_rico" : "linear",
  "inversion_height_rico" : "linear",
  "cl_acnv25_rico" : "linear",
  "cl_accr25_rico" : "linear",
}

yscaledict = {
  "thl" : "linear",
  "00rtot" : "linear",
  "rliq" : "linear",
  "prflux" : "linear",
  "cl_nc" : "linear",
  "cl_nr" : "linear",
  "clfrac" : "linear",
  "wvar" : "linear",
  "w3rd" : "linear",
  "sat_RH" : "linear",
  "rad_flx" : "linear",
  "lwp" : "linear",
  "rwp" : "linear",
  "er" : "linear",
  "wvarmax" : "linear",
  "surf_precip" : "linear",
  "acc_precip" : "linear",
  "cloud_base" : "linear",
  "gccn_rw_cl" : "linear",
  "non_gccn_rw_cl" : "linear",
  "base_prflux_vs_clhght" : "linear",
  "base_prflux_vs_clhght number of occurances" : "linear",
  "cl_gccn_conc" : "log",
  "clb_bigrain_mean_rd" : "linear", 
  "clb_bigrain_mean_kappa" :  "linear", 
  "clb_bigrain_mean_conc" :  "linear", 
  "clb_bigrain_mean_inclt" :  "linear", 
  "clb_bigrain_mean_gccn_fraction" :  "linear", 
  "cloud_cover_rico" : "linear",
  "min_cloud_base_rico" : "linear",
  "inversion_height_rico" : "linear",
  "cl_acnv25_rico" : "linear",
  "cl_accr25_rico" : "linear",
}

xlimdict_profs = {
  "thl" : None,
  "00rtot" : None,
  "rliq" : None,
  "prflux" : None,#(0,20),
  "cl_nc" : None,#(0,90),
  "cl_nr" : None,#(0,90),
  "clfrac" : None,
  "wvar" : None,
  "w3rd" : None,
  "sat_RH" : None,
  "rad_flx" : None,
  "gccn_rw_cl" : None,#(0,90),
  "non_gccn_rw_cl" : None,#(0,12),
  "base_prflux_vs_clhght" : None,#(1,10000)
  "base_prflux_vs_clhght number of occurances" : None#(1,10000)
}

ylimdict_profs = {
  "thl" : (0,3000),
  "00rtot" : (0,3000),
  "rliq" : (0,3000),
  "prflux" : (0,3000),
  "cl_nc" : (0,3000),
  "cl_nr" : (0,3000),
  "clfrac" : (0,3000),
  "wvar" : (0,3000),
  "w3rd" : (0,3000),
  "sat_RH" : (0,3000),
  "rad_flx" : (0,3000),
  "gccn_rw_cl" : (0,3000),
  "non_gccn_rw_cl" : (0,3000),
  "base_prflux_vs_clhght" : (0,2500),
  "base_prflux_vs_clhght number of occurances" : (0,2500)
}

xlimdict_series = {
  "clfrac" : (0,10),
  "cl_nc" : (0,10),
  "cl_nr" : (0,10),
  "lwp" : (0,10),
  "rwp" : (0,10),
  "er" : (0,10),
  "wvarmax" : (0,10),
  "surf_precip" : (0,10),
  "acc_precip" : (0,10),
  "cl_gccn_conc" : (0,10),
  "cloud_base" : (0,10),
  "clb_bigrain_mean_rd" : (0,10), 
  "clb_bigrain_mean_kappa" :  (0,10), 
  "clb_bigrain_mean_conc" :  (0,10), 
  "clb_bigrain_mean_inclt" :  (0,10), 
  "clb_bigrain_mean_gccn_fraction" :  (0,10),
  "cloud_cover_rico" : (0,10),
  "min_cloud_base_rico" : (0,10),
  "inversion_height_rico" : (0,10),
  "cl_acnv25_rico" : (0,10),
  "cl_accr25_rico" : (0,10),
}

ylimdict_series = {
  "clfrac" : None,
  "cl_nc" : (-1,155),
  "cl_nr" : (-0.02,.55),
  "lwp" : (-1,45),
  "rwp" : (-1,15),
  "er" : None,
  "wvarmax" : None,
  "surf_precip" : (-0.05,1.2),
  "acc_precip" : (-0.002,0.08),
  "cl_gccn_conc" : None,#(1e-6, 1),
  "cloud_base" : None,
  "clb_bigrain_mean_rd" : None, 
  "clb_bigrain_mean_kappa" :  None, 
  "clb_bigrain_mean_conc" :  None, 
  "clb_bigrain_mean_inclt" :  None, 
  "clb_bigrain_mean_gccn_fraction" :  None, 
  "cloud_cover_rico" : (-0.1,0.6),
  "min_cloud_base_rico" : (-1,750),
  "inversion_height_rico" : (1000,3000),
  "cl_acnv25_rico" : None,
  "cl_accr25_rico" : None,
}
