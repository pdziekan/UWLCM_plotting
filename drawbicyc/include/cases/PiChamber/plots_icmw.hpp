#pragma once

const std::vector<std::string> series_PiChamberICMW({
"T"
,"Qv"
,"RH"
,"S_drop"
,"LWC_gm-3"
,"IWC_gm-3"
,"N_drop"
,"N_ice"
,"N_aerosol"
,"N_drop_removal"
,"N_ice_removal"
//,"M_drop_removal"
//,"M_ice_removal"
,"disp_r"
,"r_dmean"
,"r_imean"
,"Sigma2_S"
,"Sigma2_T"
,"Sigma2_Qv"
,"Sigma2_S_drop"
,"epsilon"
,"tot_tke"
,"H_flux_t"
,"H_flux_b"
,"qv_flux_t"
,"qv_flux_b"
});

std::vector<std::string> profs_PiChamberICMW({
"T"
,"rv"
,"RH"
,"S_drop"
,"Sigma2_S"
,"Sigma2_S_drop"
,"N_c"
,"r_mean1"
}); // rtot has to be first


std::vector<std::string> fields_PiChamberICMW({
  /*
"rl", "nc"
 "rr", "nr"
//"ef", "na", 
"th", "rv",     
"u", "w", 
//"sd_conc",//, "r_dry", 
//"RH", "supersat"
//"lib_pres", "lib_temp"
*/
});
