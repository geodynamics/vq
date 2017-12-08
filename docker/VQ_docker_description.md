# Virtual Quake
Virtual Quake (VQ) is a physics-based simulation of fault stress accumulation, stress transfer, and rupture for creating synthetic seismic histories. VQ can simulate large, complex fault geometries for thousands of years on desktop-level hardware, and also includes multiprocessing capabilities.

## Initial Setup
Install Docker: <https://www.docker.com/products/docker>

Once the Docker client is running on your machine, download the Virtual Quake image with

`docker pull geodynamics/virtualquake:3.1.1`

## Running a container
First, determine a directory on your computer (The "host" computer) in which to store all simulation files. These will include any fault model files, simulation outputs, and PyVQ plots.

To allow the VQ container to access these files, include the full path of the chosen directory in the `docker run` command.  Since Docker works slightly differently on Mac and Linux, use the appropriate command for your machine:

For Linux: `docker run --rm -v <full path to host directory>:/home/virtualquake/external_volume -e HOST_UID=$(id -u) -e HOST_GID=$(id -g) -it geodynamics/virtualquake:3.1.1`
For Mac: `docker run --rm -v <full path to host directory>:/home/virtualquake/external_volume  -it geodynamics/virtualquake:3.1.1`

## Running Simulations
After executing the `docker run` command, you'll be dropped into a terminal running inside the container, in a directory with access to all the files in the host directory you chose.  The VQ fault mesher, simulator, and PyVQ analysis tools are located in the `~/vq-3.1.0` directory.

To run the mesher:
`~/vq-3.1.0/build/src/mesher <options>`

VQ simulation:
`~/vq-3.1.0/build/src/vq ./<parameter file>`

PyVQ tools:
`python ~/vq-3.1.0/PyVQ/pyvq/pyvq.py <options>`

For guidance on creating fault models, running simulations, and performing analyses (as well as available options and parameters for these tools), please reference the Virtual Quake manual, available at 
<https://geodynamics.org/cig/software/vq/>

