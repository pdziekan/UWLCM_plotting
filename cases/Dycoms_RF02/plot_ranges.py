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
  "er" : "linear",
  "wvarmax" : "linear",
  "surf_precip" : "linear",
  "acc_precip" : "linear",
  "cloud_base" : "linear",
  "gccn_rw_cl" : "linear",
  "non_gccn_rw_cl" : "linear",
  "base_prflux_vs_clhght" : "linear",
  "cl_gccn_conc" : "linear",
  "clb_bigrain_mean_rd" : "linear",
  "clb_bigrain_mean_kappa" :  "linear",
  "clb_bigrain_mean_conc" :  "linear",
  "clb_bigrain_mean_inclt" :  "linear",
  "clb_bigrain_mean_gccn_fraction" :  "linear",
  "cloud_cover_dycoms" :  "linear"
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
  "er" : "linear",
  "wvarmax" : "linear",
  "surf_precip" : "linear",
  "acc_precip" : "linear",
  "cloud_base" : "linear",
  "gccn_rw_cl" : "linear",
  "non_gccn_rw_cl" : "linear",
  "base_prflux_vs_clhght" : "log",
  "cl_gccn_conc" : "log",
  "clb_bigrain_mean_rd" : "linear",
  "clb_bigrain_mean_kappa" :  "linear",
  "clb_bigrain_mean_conc" :  "linear",
  "clb_bigrain_mean_inclt" :  "linear",
  "clb_bigrain_mean_gccn_fraction" :  "linear",
  "cloud_cover_dycoms" :  "linear"
}

xlimdict_profs = {
  "thl" : None,#(288.2,289.2),
  "00rtot" : None,#(9,10.4),
  "rliq" : (0,0.7),
  "prflux" : None,#(0,70),
  "cl_nc" : None,#(0,120),
  "cl_nr" : None,#(0,120),
  "clfrac" : None,
  "wvar" : (-0.1,0.8),
  "w3rd" : (-0.15,.15),
  "sat_RH" : (-5,1),
  "rad_flx" : None,
  "gccn_rw_cl" : None,#(0,40),
  "non_gccn_rw_cl" :None,# (0,7),
  "base_prflux_vs_clhght" : None,#(1,10000)
  "cloud_cover_dycoms" :  None,
}

ylimdict_profs = {
  "thl" : (0,1.2),
  "00rtot" : (0,1.2),
  "rliq" : (0,1.2),
  "prflux" : (0,1.2),
  "cl_nc" : (0,1.2),
  "cl_nr" : (0,1.2),
  "clfrac" : (0,1.2),
  "wvar" : (0,1.2),
  "w3rd" : (0,1.2),
  "sat_RH" : (0,1.2),
  "rad_flx" : (0,1.2),
  "gccn_rw_cl" : (0,1.2),
  "non_gccn_rw_cl" : (0,1.2),
  "base_prflux_vs_clhght" : (0,2500),
}

xlimdict_series = {
  "clfrac" : (1,6),
  "cl_nc" : (1,6),
  "cl_nr" : (1,6),
  "lwp" : (1,6),
  "er" : (1,6),
  "wvarmax" : (1,6),
  "surf_precip" : (1,6),
  "acc_precip" : (1,6),
  "cloud_base" : (1,6),
  "cl_gccn_conc" : (1,6),
  "clb_bigrain_mean_rd" : (1,6),
  "clb_bigrain_mean_kappa" :  (1,6),
  "clb_bigrain_mean_conc" :  (1,6),
  "clb_bigrain_mean_inclt" :  (1,6),
  "clb_bigrain_mean_gccn_fraction" :  (1,6),
  "cloud_cover_dycoms" :  (1,6),
}

ylimdict_series = {
  "clfrac" : None,
  "cl_nc" : (0,150),
  "cl_nr" : None,
  "lwp" : (50, 175),
  "er" : (0,1.25),
  "wvarmax" : (0,0.7),
  "surf_precip" : (-0.05,1.25),#(-0.01,0.5),
  "acc_precip" : None,#(0,0.07),
  "cloud_base" : (400,650),
  "cl_gccn_conc" : None,#(1e-10,1e-0)
  "clb_bigrain_mean_rd" : None,
  "clb_bigrain_mean_kappa" :  None,
  "clb_bigrain_mean_conc" :  None,
  "clb_bigrain_mean_inclt" :  None,
  "clb_bigrain_mean_gccn_fraction" :  None,
  "cloud_cover_dycoms" :  None,
}


