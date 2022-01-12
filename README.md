Step 0: Building the libraries needed for XSTAR.

The version of XSTAR used by PTRANSX requires a handful of libraries which are
built as a matter of course in the full installation of HEASoft. The official
documentation for how to install HEASoft will get you 80% of the way there.
However, there are some key differences---therefore, we detail the process
thoroughly here.

First, go to https://heasarc.gsfc.nasa.gov/lheasoft/download.html to download
the latest version of the HEASoft source tarball. Select the radio button for
"Source Code" then check the box for PC - Linux - CentOS (the OS that nearly all
clusters use). Scrolling down, next check "General-Use FTOOLS" and "XSTAR."
After clicking Submit, a tarball with a name like heasoft-6.29src.tar.gz will be
downloaded.

Copy this tarball to a workspace on the target HPC machine. Preferably a
non-scratch storage area; I typically name the parent directory "xstar_install."
With the tarball in this directory, issue the command:

gunzip -c heasoft-6.29src.tar.gz | tar xf -

or similar, depending on the downloaded HEASoft version number. Next, navigate
to the newly-created heasoft-6.29/BUILD_DIR (or similar). Here, run the
configure script with the following important options:

./configure --disable-x --enable-readline

Then issue:

make

This will take some time. The libraries we will need later on live in

xstar_install/heasoft-6.29/heacore/BLD/x86_64-pc-linux-gnu-libc2.17/lib

and

xstar_install/heasoft-6.29/ftools/BLD/x86_64-pc-linux-gnu-libc2.17/lib

Step 1: Cloning the repo.

Quick start: just issue

git clone https://github.com/kinchb/ptransx.git

to clone the repo, and you can get started right away.

If you'd like to contribute, you'll probably want to get added to the repo on
github. Contact brooks.e.kinch@gmail.com to make that happen.

Step 2: Creating the Compton data files.

PTRANSX and Pandurata both depend on pre-made tables for the
relativistically-correct Compton scattering photon frequency redistribution
function. These are made by the codes in the compy sub-directory, which contains
its own documentation and examples. (Please see compy's README file.) These
tables are too large to package with the code, but they need only be made once
per target machine. So let's make them now.

First, copy the whole compy directory to somewhere on a fast filesystem (like a
lustre scratch workspace). Next, examine the build script, build.sh. The
specific modules loaded are for MARCC, and these might need to be changed for a
different target machine. Compy has few dependencies: you need a C compiler (in
the example, we use Intel's), an MPI implementation (e.g., Intel's), and a
Python3 interpreter. Furthermore, we need relatively stock Python modules:
numpy, scipy, h5py, and mpi4py. The example build script loads the right SLURM
modules and installs the needed Python modules before building (in place) the
compy module. Note that building mpi4py requires an extant MPICC environment
variable set to the path of the MPI-enabled C compiler (wherever mpicc points).
This is also set in the build script. For most users, all that should need to be
changed is the "module load" commands, as indicated.

In the copied compy directory, run:

./build.sh

If successful, this creates a file like
"_mc_compy.cpython-36m-x86_64-linux-gnu.so" and a "build" directory. If we run a
Python script from within the compy directory (or add it to Python's
module-searching path), we can access compy utilities like we would any other
Python modules, with an import "import _mc_compy" command.

The "master" Compton data table is created by the p_mc_comp.py script, which can
be run in parallel by the make_table.sh script. Be sure to modify make_table.sh
with directives appropriate to your HPC platform. The example script provided is
for MARCC. Once you've modified the script appropriately, run:

sbatch make_table.sh

If successful, you will have a (large, ~18G) file in this directory called
"compton_data.h5." This file is accessed directly by PTRANSX, but Pandurata
requires a smaller version generated from it by the "pandurata_prep.py" script.
In the same directory, run:

python3 pandurata_prep.py compton_data.h5

to generate (the much smaller files) "compton_data_pan.h5" and
"compton_data_hr.h5." I suggest copying all three files to somewhere on
non-purged disk storage, just in case---but because they are accessed often in
the course of post-processing, they should live also on the fastest disk space
available.

Note that generating compton_data.h5 can take a very long time. The example
script allocates 20 nodes for 24 hours. If this is insufficient, p_mc_comp.py is
smart enough to pick up where it left off. Just re-issue the sbatch command with
the partially-complete compton_data.h5 in the directory until the table is
complete. On MARCC, this took about 3 days.

Step 3: Creating a "run" directory.

In the base PTRANSX repo directory, there is a JSON file that contains all of
the machine-specific information both to build and run PTRANSX+Pandurata. The
example case given is marcc.json. Perusing it, we see that it contains entries
for, e.g., the directories for the XSTAR libraries described in step 0, the
directories where the Compton data tables described in step 1 live (on scratch,
ideally), the modules which need to be loaded (and in the right order), and so
on. The idea is that *everything* machine-specific is encapsulated here---and
*only* here. We might have a frontera.json at some point, which would have an
identical structure but, of course, different entries.

The create_run.py script is used to create a PTRANSX run directory. Each run
directory performs the complete post-processing of a single HARM3D snapshot. I
recommend making a symlink to this script in the directory in which you want the
run directory to live. For example,

cd where_all_my_ptransx_runs_will_live

ln -s path_to_ptransx_repo/ptransx/create_run.py ./

Then, create a single run directory like so:

python3 create_run.py <config_file.json> <run_dir_name> <harm3d_snapshot.h5>

For example:

python3 create_run.py ./ptransx/marcc.json test_run test_snapshot.h5

The script will run even with no snapshot file provided, so we can test that the
machine-specific JSON file yields a working environemnt. Navigate to the
just-created directory. Notice that a symlink config.json (-> marcc.json) has
been made. The codes expects a "config.json" in the base run directory at build
and runtime. Now run the test script test.sh, like so:

./test.sh

If you get the following:

"Congratulations! You can now run PTRANSX+Pandurata!"

Then you are good to go. Otherwise, a hopefully informative error message will
print and we can debug the machine-specific JSON file together.
