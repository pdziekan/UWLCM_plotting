**Cloud Simulations with Particle-based Microphysics: Mean and Variance of Precipitation in the AON Collision-Coalescence Algorithm**

## **Running simulations**

To run any simulation please use the following line
	~/UWLCM/build_pd/uwlcm --outfreq=240 --nt=21600 --dt=0.5 --nx=121 --ny=0 --nz=101 --case=cumulus_congestus --micro=lgrngn --outdir=~/PROVIDE A NAME OF FILE_out_lgrngn --backend=multi_CUDA --sgs=0 --sd_conc_large_tail=1 --turbo_adve=0 --sd_conc={PROVIDE N_SD} --cond=1 --sstp_cond=5 --coal=1 --sstp_coal=5 --rng_seed_init=999 --rng_seed=2 --piggy=1 --vel_in=~/PATH TO/velocity_out.dat


To run any simulation one have to build [libcloudphxx](https://github.com/igfuw/libcloudphxx), [libmpdataxx](https://github.com/igfuw/libmpdataxx) and [UWLCM](https://github.com/igfuw/UWLCM). We suggest to compile and run simulation within singaularity image image_nowe.sif. 

## **Plotting results**

To obtain results as in paper one have to creat a time series of obtained results. To do so one have to build [UWLCM_plotting](https://github.com/igfuw/UWLCM_plotting), also in the singularity image. 

To run time series calculation use the following line
	~/UWLCM_plotting/drawbicyc/build/drawbicyc --dir=~/PROVIDE A NAME OF FILE_ --micro=lgrngn --type=cumulus_congestus --fields=0 --series=1 --field_plotfreq=240 --prof_start=120 --prof_end=10800

Then please go to ~/UWLCM_plotting/papers/AON_CC and use Python programs to reakreat figures from paper. Please provide PATH to fiels and PATH where to save Figs.
	pytohn3 Fig_6.py
	python3 prep_for_fig_7_8_A1.py
	python3 Fig_7.py
	python3 Fig_8.py
	python3 Fig_A1.py
