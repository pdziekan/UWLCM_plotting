#pragma once

const std::vector<std::string> series_dycoms({
// "clb_bigrain_mean_rd",
// "clb_bigrain_mean_kappa",
// "clb_bigrain_mean_conc",
// "clb_bigrain_mean_inclt",
//, "clb_bigrain_mean_gccn_fraction"

 "cl_nr",
"cl_acnv25_dycoms",
"cl_accr25_dycoms",
"wvarmax", "cloud_cover_dycoms", "lwp", "er",
 "surf_precip",
////"mass_dry",
 "acc_precip",
 "cl_nc",
 "cloud_base_dycoms"
// ,"cl_gccn_meanr"
// ,"cl_avg_cloud_rad"
// ,"cl_gccn_conc", "gccn_conc"
// "cl_rd_geq_0.8um_conc", "rd_geq_0.8um_conc"
// ,"cl_non_gccn_conc", "non_gccn_conc", "cl_gccn_to_non_gccn_conc_ratio"
// "cl_sd_conc", "cl_sd_conc_std_dev",
// "tot_water"
});

std::vector<std::string> profs_dycoms({
"00rtot",
 "rliq",
//"thl", "wvar", 
//"w3rd",
 "prflux"
//,"clfrac"
//, "N_c", 
//,"cl_nc"
//,"sat_RH"
//,"rad_flx"
//, "actrw_rw_cl"
//,"rd_geq_0.8um_conc"
//,"rd_lt_0.8um_conc"
//,"sd_conc"
//, "non_gccn_rw_cl"
//, "gccn_rw_cl_down"
//, "non_gccn_rw_cl_down"
//, "nc_up" 
//,"sat_RH_up"
//, "act_conc_up" 
//, "nc_down" 
}); // rtot has to be first

std::vector<std::string> fields_dycoms({
"rl", "nc",
"rd_geq_0.8um_conc",
"rd_lt_0.8um_conc"
// "rr", "nr",
//"ef", "na", 
//"th", "rv",     
//"u", "w", 
//"sd_conc",//, "r_dry", 
//"RH", "supersat",
//"lib_pres", "lib_temp"
//"gccn_conc",
//"gccn_mean_rw"

});
