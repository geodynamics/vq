# -*- coding: utf-8 -*-
import sys
import os.path
import readline
readline.parse_and_bind("tab: complete")

print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
print "\nThis script will help you generate a Virtual Quake parameter file.\n"
print "You will be asked about commonly modified parameters.\n"
print "If you would like to make further customizations, please consult the user manual and the example parameter file to make manual adjustements to your paramter file.\n"
print "\n"
print "After running this script, move the \"generated_parameters.d\" file to the location of your VQ run.\n"
print "A simulation can then be initiated with $[Path to the vq program] ./generated_parameters.d"
print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
print "\n"

SIM_LOCATION = raw_input("What is the full path of the full path of the directory where you'll be running this simulation?\n")
while not os.path.isdir(SIM_LOCATION):
        SIM_LOCATION = raw_input("Can't find directory, please input the full path of the directory where you'll be running this simulation\n")

year_input = raw_input("\nHow many years would you like to simulate?\n")
while True:
    try:
        END_YEAR = float(year_input)
    except:
        year_input = raw_input("Please enter a number\n")
        continue
    else:
        break

greens_exist = raw_input("\nHave you already computed the Green's functions? (y/n)\n")
while greens_exist not in ['y','Y','n','N']:
    greens_exist = raw_input("Please enter \'y\' or \'n\'\n")
if greens_exist in ['y', 'Y']:
    GREENS_METHOD = "file"
    GREENS_IOMETH = "sim.greens.input"
    GREENS_PATH = raw_input("\nWhat is the path (relative to your simulation directory) of the file where the Green's functions are stored?\n")
    while not os.path.isfile(SIM_LOCATION+GREENS_PATH):
        GREENS_PATH = raw_input("Can't find file, please input the path (relative to your simulation directory) of the Green's function file\n")
else:
    GREENS_METHOD = "standard"
    GREENS_IOMETH = "sim.greens.output"
    GREENS_PATH = raw_input("\nPlease enter the path (relative to your simulation directory) of the file where you'd like to save the Green's functions\n")
    while not os.path.isdir(os.path.dirname(SIM_LOCATION+GREENS_PATH)):
        GREENS_PATH = raw_input("Can't find directory, please input the path (relative to your simulation directory) where you'd like to save the Green's functions\n")

MODEL_FILE = raw_input("\nWhat is the path (relative to your simulation directory) of the fault model file for this simulation?\n")
while not os.path.isfile(SIM_LOCATION+MODEL_FILE):
        MODEL_FILE = raw_input("Can't find file, please input the path (relative to your simulation directory) of the model file\n")

if '.txt' in MODEL_FILE:
    MODEL_TYPE = 'text'
elif '.h5' in MODEL_FILE:
    MODEL_TYPE = 'hdf5'
else:
    MODEL_TYPE = raw_input("\nWhat is the file type for your fault model file? (hdf5 or text)\n")
    while MODEL_TYPE not in ['hdf5', 'text']:
        MODEL_TYPE = raw_input("Please enter \'hdf5\' or \'text\'\n")

EVENTS_FILE = raw_input("\nWhat is the path (relative to your simulation directory) of the file where you would like to save the events for this simulation?\n")
while not os.path.isdir(os.path.dirname(SIM_LOCATION+EVENTS_FILE)):
    EVENTS_FILE = raw_input("Can't find directory, please input the path (relative to your simulation directory) where you'd like to save the events\n")



parameter_list = ['sim.time.end_year', 'sim.system.progress_period', 'sim.greens.use_normal', 'sim.greens.method', GREENS_IOMETH, 'sim.file.input',
                  'sim.file.input_type', 'sim.file.output_event', 'sim.file.output_event_type', 'sim.friction.dynamic', 'sim.friction.dynamic_stress_drops',
                  'sim.greens.shear_offdiag_min', 'sim.greens.shear_offdiag_max', 'sim.greens.normal_offdiag_min', 'sim.greens.normal_offdiag_max']

value_list = [END_YEAR, 5, 'true', GREENS_METHOD, GREENS_PATH, MODEL_FILE, MODEL_TYPE, EVENTS_FILE, 'hdf5', 0.2, 'true',
              -745348, 3590676, -197265.50322689, 202035.55337244]


param_file = open("generated_parameters.d", "w")

print "Parameters saved as \"generated_parameters.d\".  Remember to move this file to your simulation location before running"


for i, param in enumerate(parameter_list):
    param_file.write("{:<34}= {}\n\n".format(param, value_list[i]))

param_file.close()