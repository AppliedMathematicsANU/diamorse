1) Packages required for `make` or `make main` (just the C++ programs):

    gcc-g++ (Devel)
    make (Devel)
    git (Devel)
    time (Utils)


2) Additional packages required for `make python` (Python wrappers):

    python-numpy (Python)
    python-cython (Python)

After running `make python`, rename the file `bin/MorseAnalysis.so` to `bin/MorseAnalysis.dll`.


3) Additional packages required in order to make plots in python (file output only):

    libX11-devel (Libs)
    libfreetype-devel (Libs)

With those packages installed, run Cygwin as administrator and issue the following commands:

    easy_install-2.7 pip

(replace the "2.7" with the Python version you're using.)

    pip install matplotlib

The script `plot_basins.py` in the `python/` directory has an option `-o` for writing the output to a file. The script   `plot_persistence.py` does not have such an option yet, as it is essentially just a template for writing custom scripts. The only change needed is to replace the line `plt.show()` at the end of the file with something like `plot.savefig("figure.png")`.
