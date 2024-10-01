\page Example_003_dld Example #3: Configuring particles on a streakline to be injected into a microfluidic DLD device
\tableofcontents 

# Overview

The files for this example can be found in `pt/docs/examples/dld_03_500nm`.
This example considers a microfluidic DLD structure with a periodicity of $Np=8$ and a circular pillar array of diameter of $1~\mu m$ and pitch of $2~\mu m$.

The continuous phase solution fields have already been obtained through solving the fluid flow equations of motion by using Ansys Fluent package. 

# Hierarchy of directories and files
- The `./inp` directory is used to store the mesh and solution files. Make sure you copy the mesh and solution files from `pt/docs/examples/dld_01/inp` into this directory.
- The `./particle` directory is created before running the simulation and will be used by the library to store the particle trajectory files.
- The `./report` directory is created before running the simulation and will be used by the library to store the report files as needed.

# Configuration file

Navigate to the project directory. 
Inspect the configuration file named `config` and review the comments included in the file that aims at explaining the structure of the file and the purpose of using each parameter.

The configuration file is set up for the case of injecting 321 particles of diameter $500~nm$ and density of $1,000~kg/m^3$ on a streakline connecting two points of $(1~\mu m,~0~\mu m)$ and $(1~\mu m,~32~\mu m)$, resulting in a $100~nm$ distance between initial position of particles.

# Running particle tracking simulation

While you are in the aformentioned directory in your terminal, run:

```
pt -f config
```

The particle tracking will be completed in about a minute or so based on your computer specifications.
The particle trajectory files shoule be accessible from the `./particle` directory.
Here are some of the results:

> [!NOTE]
> Feel free to look at similar simulations for other particle diameters that can be found in `pt/docs/examples/dld_03_300nm`, `pt/docs/examples/dld_03_400nm`, and `pt/docs/examples/dld_03_600nm`. In particular, it can be discerned how particle motion changes from zig-zag maode of small particles to bumping mode of large particles. 

<img src="Example_003_dld_res_00.png" alt="Image" style="width:75%;"/>

\image latex Example_003_dld_res_00.png "Trajectory of particles obtained by using the library." width=\textwidth