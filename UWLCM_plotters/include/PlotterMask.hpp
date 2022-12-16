#pragma once
#include "PlotterMicro.hpp"

enum class mask_type_t{Rico11, Dycoms_rf02};

template<int NDims>
class PlotterMask : public PlotterMicro<NDims> 
{
  public:
  using parent_t = PlotterMicro<NDims>;
  using arr_t = typename parent_t::arr_t;

  private:

  mask_type_t mask_type;
  arr_t mask;

  void calc_mask(int at)
  {
    if(mask_type == mask_type_t::Rico11)
    {
      mask = this->load_rc_timestep(at);
      mask = iscloudy_rc_rico(mask);
    }
    else if(mask_type == mask_type_t::Dycoms_rf02)
    {
      mask = this->load_nc_timestep(at);
      mask = iscloudy_nc_dycoms(mask);
    }
    else throw std::runtime_error("Invalid mask type");
    mask(this->hrzntl_slice(0)) = 0; // cheat to avoid occasional "cloudy" cell at ground level due to activation from surf flux
  }

  //lgrngn_droplet_prefix

  private:
  // helper function that calculates staistics (mean and std_dev) of a field only in cloudy cells
  std::pair<double, double> cloud_hlpr(arr_t arr, int at)
  {
    std::pair<double, double> res;

    // apply mask
    calc_mask(at);
    arr *= mask; 
    
    if(blitz::sum(mask) > 0.) 
      res.first = blitz::sum(arr) / blitz::sum(mask); 
    else
      res.first = 0.; 

    arr = pow(arr - res.first, 2); 
    arr *= mask; // apply filter
    if(res.first>0)
      res.second = sqrt(blitz::sum(arr) / blitz::sum(mask)); 
    else
      res.second = 0.;

    return res;
  }

  public:


  arr_t get_mask(int at)
  {
    calc_mask(at);
    return mask;
  }

  // functions for diagnosing statistics
  // mean and std dev [g/kg] of the mixing ratio of activated dropelts in cloudy cells (characteristics of the spatial distribution at this timestep)
  std::pair<double, double> cloud_ract_stats_timestep(int at)
  {
    // read activated droplets mixing ratio 
    arr_t ract(this->load_ract_timestep(at));
    ract *= 1e3; // turn it into g/kg
    return cloud_hlpr(ract, at);
  }

  // mean and std_dev of concentration of activated droplets in cloudy cells [1/cm^3] (characteristics of the spatial distribution at this timestep)
  std::pair<double, double> cloud_conc_stats_timestep_hlpr(int at, std::string lgrngn_prefix, std::string blk_2m_1, std::string blk_2m_2 = "")
  {   
    arr_t conc;
    // read concentration of activated droplets
    if(this->micro == "blk_1m") return {0,0};
    // TODO: fix stupid copying of arrays
    else if(this->micro == "lgrngn") {
      arr_t tmp(this->h5load_timestep(lgrngn_prefix+"_mom0", at));
      conc.resize(tmp.shape());
      conc = tmp;
    }
    else if(this->micro == "blk_2m") {
      arr_t tmp(arr_t(this->h5load_timestep(blk_2m_1, at)));
      if(blk_2m_2 != "")
        tmp += arr_t(this->h5load_timestep(blk_2m_2, at));
      conc.resize(tmp.shape());
      conc = tmp;
    }
    this->multiply_by_rhod(conc); // b4 it was specific moment
    conc /= 1e6; // per cm^3
    return cloud_hlpr(conc, at);
  } 

  // mean and std_dev of concentration of activated droplets in cloudy cells [1/cm^3] (characteristics of the spatial distribution at this timestep)
  std::pair<double, double> cloud_actconc_stats_timestep(int at)
  {   
    return cloud_conc_stats_timestep_hlpr(at, "act_rw", "nc", "nr");
  } 

  std::pair<double, double> cloud_cloudconc_stats_timestep(int at)
  {   
    return cloud_conc_stats_timestep_hlpr(at, "cloud_rw", "nc");
  } 

  std::pair<double, double> cloud_rainconc_stats_timestep(int at)
  {   
    return cloud_conc_stats_timestep_hlpr(at, "rain_rw", "nr");
  } 


  // mean and std_dev of number of SDs in cloudy cells (characteristics of the spatial distribution at this timestep)
  std::pair<double, double> cloud_sdconc_stats_timestep(int at)
  {   
    if(this->micro == "blk_1m") return {0,0};
    arr_t sdconc(this->h5load_timestep("sd_conc", at));
    return cloud_hlpr(sdconc, at);
  }

  // mean and std_dev of number of activated SDs in cloudy cells (characteristics of the spatial distribution at this timestep)
  std::pair<double, double> cloud_sdconc_act_stats_timestep(int at)
  {   
    if(this->micro == "blk_1m") return {0,0};
    arr_t sdconc_act(this->h5load_timestep("sd_conc_act", at));
    return cloud_hlpr(sdconc_act, at);
  }

  // mean and std_dev of mean radius of activated droplets in cloudy cells [um] (characteristics of the spatial distribution at this timestep)
  std::pair<double, double> cloud_meanr_stats_timestep(int at)
  {   
    if(this->micro == "blk_1m") return {0,0};
    // read act drop 0th raw moment / mass [1/kg]
    arr_t act0th(this->h5load_timestep("actrw_rw_mom0", at)); 
    // read act drop 1st raw moment / mass [um/kg]
    arr_t act1st(this->h5load_timestep("actrw_rw_mom1", at) * 1e6);
    // calculate mean radius, store in act1st
    act1st = where(act0th > 0, act1st / act0th, 0.);
    return cloud_hlpr(act1st, at);
  }

  // mean and std_dev of std_dev of radius of activated droplets in cloudy cells [um] (characteristics of the spatial distribution at this timestep)
  std::pair<double, double> cloud_stddevr_stats_timestep(int at)
  {   
    if(this->micro == "blk_1m") return {0,0};
    // read act drop 0th raw moment / mass [1/kg]
    arr_t act0th(this->h5load_timestep("actrw_rw_mom0", at)); 
    // read act drop 1st raw moment / mass [um/kg]
    arr_t act1st(this->h5load_timestep("actrw_rw_mom1", at) * 1e6);
    // read act drop 2nd raw moment / mass [um^2/kg]
    arr_t act2nd(this->h5load_timestep("actrw_rw_mom2", at) * 1e12);
    // calculate stddev of radius, store in act1st
    act1st = where(act0th > 0, 
      act2nd / act0th - act1st / act0th * act1st / act0th, 0.);
    // might be slightly negative due to numerical errors
    act1st = where(act1st < 0, 0, act1st);
    act1st = sqrt(act1st);
    return cloud_hlpr(act1st, at);
  }

  
  // height [m] of the center of mass of activated droplets
  double act_com_z_timestep(int at)
  {
    arr_t ract(this->load_ract_timestep(at));
    arr_t weighted(ract.copy());
    weighted = weighted * this->LastIndex * this->map["dz"];
    if(blitz::sum(ract) > 1e-3)
      return blitz::sum(weighted) / blitz::sum(ract);
    else
      return 0.; 
  }


  //ctor
  PlotterMask(const string &file, const string &micro, const mask_type_t _mt):
    parent_t(file, micro),
    mask(this->tmp_ref.shape()),
    mask_type(_mt)
  {
  }
};

