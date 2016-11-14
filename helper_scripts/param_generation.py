# -*- coding: utf-8 -*-
import sys
import os.path
import readline
readline.parse_and_bind("tab: complete")

print "\nThis script will help you generate a Virtual Quake parameter file.\n"
print "You will be asked about commonly modified parameters.\n"
print "If you would like to make further customizations, please consult the user manual and the example parameter file to make manual adjustements to your paramter file.\n"
print "\n"

year_input = raw_input("How many years would you like to simulate?\n")
while True:
    try:
        END_YEAR = float(year_input)
    except:
        year_input = raw_input("Please enter a number\n")
        continue
    else:
        break

greens_exist = raw_input("Have you already computed the Green's functions? (y/n)\n")
while greens_exist not in ['y','Y','n','N']:
    greens_exist = raw_input("Please enter \'y\' or \'n\'\n")
if greens_exist in ['y', 'Y']:
    GREENS_METHOD = "file"
    GREENS_PATH = raw_input("What is the full path of the file where the Green's functions are stored?\n")
    while not os.path.isfile(GREENS_INPUT):
        GREENS_PATH = raw_input("Can't find file, please input the full path of the Green's function file\n")
else:
    GREENS_METHOD = "standard"
    GREENS_PATH = raw_input("Please enter the full path of the file where you'd like to save the Green's functions\n")
    while not os.path.isdir(os.path.dirname(GREENS_PATH)):
        GREENS_PATH = raw_input("Can't find directory, please input the full path where you'd like to save the Green's functions\n")

MODEL_FILE = raw_input("What is the full path of the fault model file for this simulation?\n")
while not os.path.isfile(MODEL_FILE):
        MODEL_FILE = raw_input("Can't find file, please input the full path of the model file\n")


EVENTS_FILE = raw_input("What is the full path of the file where you would like to save the events for this simulation?\n")
while not os.path.isdir(os.path.dirname(EVENTS_FILE)):
    EVENTS_FILE = raw_input("Can't find directory, please input the full path where you'd like to save the events\n")



parameter_list = ['sim.time.end_year', 'sim.system.progress_period', 'sim.greens.use_normal', 'sim.greens.method', 'sim.greens.input', 'sim.file.input',
                  'sim.file.input_type', 'sim.file.output_event', 'sim.file.output_event_type', 'sim.friction.dynamic', 'sim.friction.dynamic_stress_drops',
                  'sim.greens.shear_offdiag_min', 'sim.greens.shear_offdiag_max', 'sim.greens.normal_offdiag_min', 'sim.greens.normal_offdiag_max']

value_list = [END_YEAR, 5, 'true', GREENS_METHOD, GREENS_PATH, MODEL_FILE, 'hdf5', EVENTS_FILE, 'hdf5', 0.2, 'true',
              -745348, 3590676, -197265.50322689, 202035.55337244]


param_file = open("generated_parameters.d", "w")

print "\"Parameters saved as generated_parameters.d.\""


for i, param in enumerate(parameter_list):
    param_file.write("{:<34}= {}\n\n".format(param, value_list[i]))

param_file.close()