# Calling pyvq from the command line (all commands are one line)
 ----------------- Example usage --------------------------

 Frequency magnitude plot for all events in the sim file
 python path/to/vq/pyvq/pyvq.py --event_file path/to/sim_file.h5 
            --event_file_type 'hdf5' --plot_freq_mag 1

 Frequency magnitude plot including a b=1 line, only using EQs that produce slip
   on sections 1,2,3 (must also specify a model file)
 python path/to/vq/pyvq/pyvq.py --event_file path/to/sim_file.h5 
            --event_file_type 'hdf5' --plot_freq_mag 1 
            --model_file path/to/model.txt --use_sections 1 2 3

 Magnitude vs Rupture Area plot for all events in the sim file
 python path/to/vq/pyvq/pyvq.py --event_file path/to/sim_file.h5 
            --event_file_type 'hdf5' --plot_mag_rupt_area

 Plot prob vs time of 3.0 < Mag < 6.0 earthquakes for all events in the sim file
 python path/to/vq/pyvq/pyvq.py --event_file path/to/sim_file.h5 
            --event_file_type 'hdf5' --plot_prob_vs_t --min_magnitude 3.0 
            --max_magnitude 6.0

 Plot conditional prob vs time of Mag > 5.5 earthquakes, and plotting a weibull
   curve with beta=1.2 and tau=19
 python path/to/vq/pyvq/pyvq.py --event_file path/to/sim_file.h5 
            --event_file_type 'hdf5' --plot_cond_prob_vs_t --min_magnitude 5.5 
            --beta 1.2 --tau 19

 Plot waiting times vs time for Mag > 6.0 earthquakes
 python path/to/vq/pyvq/pyvq.py --event_file path/to/sim_file.h5 
            --event_file_type 'hdf5' --plot_waiting_times --min_magnitude 6.0
            
 Plot gravity changes for 5m slip on all elements in fault model
 python path/to/vq/pyvq/pyvq.py --model_file path/to/faults.txt --field_plot
            --field_type "gravity"


 ----- Currently working plot options (with examples where applicable) ---------
 --plot_freq_mag 1 (normal) or 2 (add a b=1 line) or 3 (add UCERF2 observed 
                           rates with errorbars) or 4 (add b=1 and UCERF2 rates)
 --plot_mag_rupt_area
 --plot_mag_mean_slip
 --plot_prob_vs_t
 --plot_prob_vs_t_fixed_dt
 --plot_cond_prob_vs_t
 --plot_waiting_times
 --min_magnitude 5.0
 --max_magnitude 7.9
 --use_sections 1 2 3 4
 --field_plot 
 --field_type 'gravity' (options: gravity, dilat_gravity, displacement, insar, 
                                  potential, geoid)
 --colorbar_max 200 (there are default values for each field_type, but e.g. to 
                     change colorbar range on gravity plot to 
                     -200 microgal < colors < 200 microgal)

