# Particle Simulation in Headless Mode

In this project we aim to port various particle based simulation
codes, for fluids such as SPH (smoothed particle hydrodynamics) or for
regular particles so that they run without graphics, dump their output
to a position file which we plot in a seperate script. So there will
be different solvers with their visual parts taken out and only the
skeletal location computation left. All particle locations will be
dumped to a CSV file, which will be processed by a Python Mayavi based
visualizer.

Dependencies

Each subproject will list its own dependencies which, as I mentioned
above, will be minimal. Compiling will be easy, simply run `make`. On
Ubuntu `apt-get install` on

```
build-essentials
```

package. For the Python visualizers the dependencies are

```
numpy
mayavi
pandas
```

Running

Unless they have own file output mechanisms, subprojects will dump its
particle locations under `tmp` in files `simsph-x.csv`,
`simsph-y.csv`, `simsph-z.csv` for each coordinate, every line will be
a time step. The `;` seperated CSV file in each line will carry the
particle location for that axis. The Python visualizer processes these
files with `viewsph.py`. This code dumps its image files, for each
time step under `/tmp/simsph`. From individual image files it is easy
to create a GIF animation using ImageMagick,

```
convert -delay 30 /tmp/simsph/*.png $HOME/Downloads/simsph.gif
```

# Sample Output

This output comes from `fluid-engine-dev` subproject, it uses over
200K particles, offline computed on CPU.

<img width="80%" src="https://drive.google.com/uc?export=view&id=1Idz5In7YwHolymWd1DVXMKMXeBhVRX4Q"/>

