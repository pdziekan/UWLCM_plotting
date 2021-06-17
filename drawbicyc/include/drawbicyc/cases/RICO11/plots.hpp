#pragma once

const std::vector<std::string> series_rico({
// "clb_bigrain_mean_rd",
// "clb_bigrain_mean_kappa",
// "clb_bigrain_mean_conc",
// "clb_bigrain_mean_inclt",
//, "clb_bigrain_mean_gccn_fraction"

"cl_acnv25_rico",
"cl_accr25_rico",
"cloud_cover_rico",
"min_cloud_base_rico",
"inversion_height_rico",
"lwp",
"rwp",
"surf_precip",
"acc_precip",
"cl_nc",
"cl_nr",
"sd_conc"

/*
 "cloud_base",
 "surf_flux_latent",
 "surf_flux_sensible"
 ,"cl_sd_conc"
//"mass_dry",
 ,"cl_gccn_conc", "gccn_conc"
 ,"cl_non_gccn_conc", "non_gccn_conc", "cl_gccn_to_non_gccn_conc_ratio"
//, "cl_gccn_meanr"
//,"cl_avg_cloud_rad"
// "cl_sd_conc_std_dev",
// "tot_water"
*/
});

std::vector<std::string> profs_rico({
"00rtot"
//,"rd_lt_0.8um_conc"
//,"rd_geq_0.8um_conc"
//,"sd_conc"
, "rliq"
//, "thl", "wvar", 
 ,"prflux"
//,"clfrac"
//,"sd_conc"
//,"cl_nc"
//,"cl_nc_up"
//,"w"
//,"u", "v"
//,"base_prflux_vs_clhght"
//, "non_gccn_rw_cl"
//, "gccn_rw_cl"
//,"sat_RH_up"
//, "N_c", 
//,"vel_div"
//, "nc_up" 
//,"sat_RH_up"
//, "act_conc_up" 
//, "nc_down" 
}); // rtot has to be first


std::vector<std::string> fields_rico({
"rl", "nc",
// "rr", "nr",
//"ef", "na", 
//"th", "rv",     
//"u", "w", 
//"sd_conc",//, "r_dry", 
//"RH", "supersat",
//"lib_pres", "lib_temp"
"rd_geq_0.8um_conc",
"rd_lt_0.8um_conc"
});
