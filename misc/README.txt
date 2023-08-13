ctrl+Z
pkill -9 commander3

PROJECT. INSTRUCTIONS/HELP:
    QUICK SUMMARY - Commander project: [Hans Kristian Eriksen]
        Here's a quick summary on how to get started on the Commander project:
        1) Choose a dataset to work with, identify the Healpix map on Lambda or elsewhere. Plot that map, and make sure it looks like a diffuse radio map. Does it have an extended Galaxy? If not, this is probably not an experiment you want to work with.
        2) Read the paper that describes the data set. Don't worry if you don't understand everything, but it's important to have a rough idea about how the data were taken.
        3) Create a noise RMS map. For many experiments, this will simply be a constant number across the full sky with typical sensitivity per pixel. If so, create a Healpix map with the same Nside as the main map, and fill the entire map with that number. Other experiments may provide the map for you.
        4) Create a beam file. This will typically be given as a FWHM in arcmin by the original paper. If so, create a tabulated beam with Healpy, and store that as a Healpix beam fits file.
        5) Create a bandpass file. This will often be given as a top hat function in the original paper. If so, create an ASCII file with two lines: (nu_low, 1) for the low boundary on the first, and (nu_high, 1) for the high boundary on the second. If you only have the band center, then choose 'delta' as your bandpass type
        6) Identify which units the map is given in.
        7) Make a mask that consists of 1 for accepted pixels and 0 for rejected pixels. If you see visual crap in the map by eye, remove those pixels as well.
        8) All maps must be in Galactic coordinates.
        Store all data products in the data directory specified in  the Commander parameter file.
        Congratulations, you are now ready to create the Commander parameter file section, and run this!


    IMPORTANT DIRECTORY INFO:
        Collective Parent directory: 
            /mn/stornext/d16/cmbco/AST9240/2022
        Collective Parameter file: 
            /mn/stornext/d16/cmbco/AST9240/2022/param_AST9240.txt

        Please create your own working directory within the parent directory, copy over the parameter file, and try to run for an iteration!

        [GIT] Also, make sure we're all on the Commander branch: 
            https://github.com/Cosmoglobe/Commander/tree/AST9240_2022

        When adding your data to the data directory, please make your own directory (e.g. simran) to store your maps in. 
            Then in your default file you create, just add that directory name in front of your filename:
                BAND_MAPFILE&&& = simran/filename.fits


    OWL INFO:
        The following owls are generally off limits:
            Owl36
            owl37 (WMAP+LFI TOD run)
            owl31 (hosts ganglia.owl.no)
            owl24 (dedicated for COMAP processing)
        
        There are issues with the lower OWLS (19 to 24), so do not install Commander on those.

    COLLECTIVE RUN:
        [Wang Wang]
            For creating our collective run: 
                Defaults directory: 
                    /mn/stornext/d16/cmbco/AST9240/2022/defaults
                Collective parameter file: 
                    /mn/stornext/d16/cmbco/AST9240/2022/param_AST9240_tmp.txt
                Common dipole file: 
                    /mn/stornext/d16/cmbco/AST9240/2022/data/md_hke_tmp.dat
                Common instrument file:
                    /mn/stornext/d16/cmbco/AST9240/2022/data/instrument_ast9240.dat

                (1) copy your defaults file to the  Defaults directory,
                (2) modify the common file to include your dipole, gain, bandpass
                (3) modify Collective parameter file to add your map, (in "our map sets" section)
        
        [Daniel Herman]
            Here is the common instrument file for adjusting your gains: 
                /mn/stornext/d16/cmbco/AST9240/2022/instrument_ast9240.dat

        [Lukas Hergt]
            (QUESTIONS)
            What would be the header for the md_hke_tmp.dat file? In other words, what do the various columns correspond to?
            Follow-up question: 
                How do you determine the values that you put in that file for your corresponding dataset?
        
            (ANSWER)
            Answer to my first question:
                monopole dipole_x dipole_y dipole_z prior_mean_mp prior_mean_dp_x prior_mean_dp_y prior_mean_dp_z std_mp std_dp

            (QUESTION)
            Does the md_hke_tmp.dat file assume a specific unit? Which one?
        
        [Daniel Herman]
            (ANSWER)
            It assumes the same unit as the corresponding band


    PAPER WRITING:
        [Daniel Herman]
            In the parent directory /mn/stornext/d16/cmbco/AST9240/2022 there is references.bib file. Please add a paper reference for your data set, preferably with a descriptive name!

        OVERLEAF: [Daniel Herman]
            Moving forward, here's a link to the overleaf for the proposed low-frequency paper: 
                https://www.overleaf.com/2234218746ftmwvfxzhpgr
            Your task is to add your data set to the data section, and write a paragraph (~4 lines) giving an overview of your experiment.




COMMANDER. CODE INSTRUCTIONS/HELP:
    IMPORTANT DIRECTORY INFO:
        [me] Commander installation dir and subdirs: 
            /mn/stornext/u3/jibran/Commander/
            /mn/stornext/u3/jibran/Commander/build_owl2528_oneapi
        [me] Bash file for RUNNING COMMANDER: 
            /mn/stornext/d16/cmbco/AST9240/2022/jibran/run_Commander.sh
        [me] VERY important data file (but what does it contain?):
            /mn/stornext/d16/cmbco/AST9240/2022/data/md_jibran.dat
        [me] Defaults lwa1 file for commander: 
            /mn/stornext/d16/cmbco/AST9240/2022/jibran/lwa1_80mhz_jibran.defaults
        [me] Params lwa1 file for commander:
            /mn/stornext/d16/cmbco/AST9240/2022/jibran/param_jibran_lwa1_80mhz.txt
        [me] Input lwa1 fit files for commander:
            /mn/stornext/d16/cmbco/AST9240/2022/data/jibran/
                    The following dir is not to be used as input dir. This is just for me to store fit files if I want: /mn/stornext/d16/cmbco/AST9240/2022/jibran/commander-files
        [me] Output data from commander:
            /mn/stornext/d16/cmbco/AST9240/2022/jibran/commander-files
        

    INSTALLATION: [Maksym Brilenkov]
        1. Log-in into login server with your user name and the password:
                $ ssh -X <user_name>@login.astro.uio.no

        2. Log-in into one of the owls:
            $ ssh -X owl<#>.uio.no
            
            For instance:
                $ ssh -X owl28.uio.no
            You can check which owl is available here:
                http://owl.uio.no/
            Just press the “Refresh Data” button to get the latest update of the resources available.
            
        3. Choose the directory where you want to store Commander and clone it there:
                $ git clone https://github.com/Cosmoglobe/Commander.git AST9240_2022

        4. Go into that folder and checkout the AST9240_2022 branch:
                $ cd AST9240_2022 && git checkout AST9240_2022

        5. Compile commander using the install_ita.sh script:
                $ ./install_ita.sh
            
            Note: Ensure to remove anaconda from your .bashrc if it is there for some reason. You can just comment it out while commander installs and then comment it back.
            
            If, it crushed with the following error message:
                [ 55%] Performing build step for 'hdf5'
                CMake Error at /mn/stornext/u3/maksymb/commander/wmap_test/build_owl3135_oneapi/subbuilds/src/hdf5-stamp/hdf5-build-RelWithDebInfo.cmake:49 (message):
                Command failed: 2

                '/usr/bin/gmake'

                See also

                    /mn/stornext/u3/maksymb/commander/wmap_test/build_owl3135_oneapi/install/logs/hdf5-build-*.log

            try to restart the install_ita.sh  script. If that didn’t work either, then talk to me.

        6. Once you installed it, you need to update your .bashrc otherwise it will not run. Go to your $HOME  directory:
            and open .bashrc  with your text editor and add there the following lines:
                export HEALPIX=<path to the commander root directory>/build_owl<>_<compiler>/install/healpix
                export COMMANDER_PARAMS_DEFAULT="<path to the commander root directory>/commander3/parameter_files/defaults/"
            Modules you need to load to run it:
                module load intel/oneapi mpi/latest icc/latest compiler-rt/latest mkl/latest

    COMMANDER SETUP for Run:
        1. INSTALL (above)

        2. Make following edits/changes/additions to relevant files:
            (a) .bashrc
                    export HEALPIX="~/name-of-commander-install-dir/build_owlXXXX_oneapi/install/healpix"
                    export COMMANDER_PARAMS_DEFAULT="~/AST9240_2022/commander3/parameter_files/defaults/"
            (b) param_jibran_lwa1_80mhz
            (c) run_Commander.sh
            (d) lwa1_80mhz_jibran.defaults
            (e) md_jibran.dat

        3. Remember to copy latest LWA1 fits files to the collective data dir (inside 'jibran/')
    
    SAMPLING THE GAIN FOR BANDS: [Daniel Herman]
        If you want to be able to sample the gain for bands:
            Put
                SAMPLE_SPECTRAL_INDICES = .true.
            Under the synchrotron component add:
                COMP_BETA_PRIOR_GAUSS_RMS&& = 0.0
            Under the thermal dust component add:
                COMP_BETA_PRIOR_GAUSS_RMS&& = 0.0
            Under the AME component add:
                COMP_NU_P_PRIOR_GAUSS_RMS&& = 0.0

    BASH FILE FOR RUNNING COMMANDER: [Duncan Watts]
        For running commander, I have a large bash file that has things like this:
            killall -9 commander3
            export COMMANDER_PARAMS_DEFAULT=$HOME"/Commander/commander3/parameter_files/defaults/"

            build=owl3135
            pfile=param_WMAP_only.txt
            dir=chains_WMAP_220825

            mkdir -p $dir
            cp $pfile $dir/param_orig.txt
            mpiexec -env I_MPI_FABRICS shm -n $n /mn/stornext/u3/duncanwa/Commander/build_$build"_oneapi/install/bin/commander3" $pfile --OUTPUT_DIRECTORY=$dir 2>&1| tee $dir/slurm.txt

        Some of these are specific output directories 
            (For example, /mn/stornext/u3/duncanwa/Commander/build_$build"_oneapi/install/bin/commander3 should point to your own compiled binary...) 
        and use a number like 
            n=8 
        to get started.

    OTHER USEFUL STUFF for Commander:
        As a check, one can set OUTPUT_SIGNALS_PER_BAND = .true.  to be able to check the signals for each band

        KILLING COMMANDER: [DUNCAN WATTS]
            My workflow for killing commander:
                ctrl+Z
                pkill -9 commander3
            you can also see the process ID on top, and do 
                kill -9 PROCESS_NUMBER
        
        DEBUGGING:
        Floating point error
            [Hans Kristian]
            I see that you have set the dipole uncertainty in md_anshul.dat to zero (last number on your line). This cannot be zero, since you then divide by zero. Set it to a small positive value.



MISC. USEFUL:
    AST9240 CompSep Course GIT Repo:
        https://github.com/Cosmoglobe/Component-Separation-Course

    Spreadsheet with names and datasets:
        https://docs.google.com/spreadsheets/d/18rGqtv0M0UkVE5vY1vWTpCKUC82tka-_XneoemcXliA/edit?usp=sharing

    Datasets:
        LAMBDA datasets in general: 
        https://lambda.gsfc.nasa.gov/product/
        polarization surveys
        https://lambda.gsfc.nasa.gov/product/foreground/fg_pol_survey.html
        diffuse foregrounds in general:
        https://lambda.gsfc.nasa.gov/product/foreground/fg_diffuse.html

    Build your own Gibbs Sampler:
        https://github.com/Cosmoglobe/Component-Separation-Course/blob/main/build_your_own_gibbs_sampler.ipynb

    Slides for the AST9240 Course:
        https://drive.google.com/drive/folders/1HQhg0rQjHOfoZunrlvBGbAV0Ppxiuyan?usp=sharing

    Data not published in a table??:
        For those of you looking for data that isn't published in a table, this website might be useful:
        https://apps.automeris.io/wpd/

    Tutorial on rotating Healpy maps:
        https://zonca.dev/2021/03/rotate-maps-healpy.html

        [STUART] One thing to keep in mind with this is that you can get some reasonably large artefacts in the rotated maps when rotating an incomplete sky in a_lm space. Usually the worst effect occurs near the edge pixels, but depending on the mask you can get relatively large aliasing of power across the whole rotated map. Of course if you are rotating only for visual purposes then this isn’t a big deal. Rotating in pixel space is generally safer, but since there isn’t a one-to-one mapping you get a smoothing of your data that you should take into account.

