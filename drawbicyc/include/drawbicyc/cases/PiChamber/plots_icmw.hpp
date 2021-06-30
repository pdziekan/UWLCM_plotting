#pragma once

const std::vector<std::string> series_PiChamberICMW({
"T"
,"Qv"
,"RH"
,"S_drop"
//,"LWC"
,"LWC_gm-3"
,"N_drop"
,"N_aerosol"
,"N_removal"
,"disp_r"
,"r_mean1"
//,"r_mean2"
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
,"Sigma2_S"
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
