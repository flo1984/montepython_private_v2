from procoli import lkl_prof

profile = lkl_prof(
  chains_dir='../../output/NEDE_Plik_Mb_Panth_Plus_dlogw_dloga/', 
  prof_param='H0', 
)

# Set the other profile parameters either through the class above 
# or as shown below 
profile.prof_max = 72
profile.prof_min = 65.
profile.processes = 6

profile.prof_incr = 0.1 # run two separate jobs with 
                        # both a + increment and a - increment

# Settings for global best fit search 
# The below are the defaults 
profile.set_global_jump_fac([1, 0.8, 0.5, 0.2, 0.1, 0.05])
profile.set_global_temp([0.3333, 0.25, 0.2, 0.1, 0.005, 0.001])

# Check the global best fit, run if necessary, 
# record this point in the likelihood profile txt file 
profile.global_min(
   run_glob_min=True,  
   N_min_steps=4000 
) 

# Print as a check 
print("Global minimum: ")
print(profile.global_ML)



