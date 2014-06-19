RELATIVE LUMINOSITY ANALYSIS
----------------------------

Installation
------------
  0. run "install" to make directories and build symlinks; you should
     have this current directory as a subdirectory of root12fms
  1. build scaler bit reader by typing "make clean" and then "make"
  2. download scaler data, and drop it in the directory ./sca2013
     (see instructions under "detailed version")
  3. download spin pattern data, and drop it in ./spinpat



Analysis Procedure
------------------
1. Download scaler files from HPSS

  - 2013 instructions
    - code contained in "hpss" subdirectory locally; on rcas found in ~/scalers2013
      - "BuildList" uses goodruns.dat (see trgmon code) to build a list of
        scaler files to retrieve from HPSS; the variable "size" controls
        how many scaler files will be listed in each HPSS request list
        - note that a board number is needed; see usage note
        - the target is set to $HOME/scratch/sca$YEAR/...
      - run "hpss_user.pl -f [request list]" to add the list of files
        and respective targets to the HPSS request queue and wait
      - monitor status with "hpss_user.pl -w"
  - 2012 instructions
    - scaler files on HPSS for run12 were saved including a UNIX timestamp, which
      makes it difficult to generate a list of files
    - first, type "hsi" and cd to /home/starsink/raw/daq/2012
    - then type "out > FILE_LIST"; this will create a file called "FILE_LIST" in 
      your current local (RCAS) directory and pipe all output from hsi to this file
    - run12 rellum uses board 12 data, so to list all the board 12 scaler files, type
      "ls */*/*_12_*.sca" and be patient
    - when the command is done, exit HPSS and check FILE_LIST for the proper output
    - execute BuildList_2012 to build lists of files to submit to the data carousel; 
      files_to_retrieve*.lst will be created, with 150 requests per file
      - submit using "hpss_user.pl -f [list]"
      - check carousel status using "hpss_user.pl -w"


2. Read the scalers scaler bit reader reader from Zilong
   (32-bit for run13 and 24-bit for run12)

  - 32-bit scaler reader (for 2013)
    (obtained from ~zchang/2013-03-scaler/codes_scaler/MyCodes/scaler2_reader_bit.c)
    - read_scalers in scalers13 directory executes scaler2_reader_bit.exe via condor
      - reads all scaler files in /GreyDisk1/sca2013
      - outputs the read files in datfiles directory
      - explanation of scaler2_reader_bit.exe (see make file for compilation)
        - reads 32-bit scaler files obtained from hpss (downloaded to /GreyDisk1/sca2013)
        - outputs corresponding file in datfiles directory with the following columns:
          - bunch crossing (bXing) number
          - BBC[0-7] (see zilong's scaler bit definitions below)
          - ZDC[0-7] 
          - VPD[0-3]
          - total bXings
  - 24-bit scaler reader (for 2012)
    - read_scalers in scalers12 directory executes sca_read_bin.o via condor
      - reads all scaler files in /GreyDisk1/sca2012
      - outputs the read files in datfiles directory
      - datfile columns
        - bunch crossing (bXing) number
        - BBC[0-7] (see zilong's scaler bit definitions below)
        - ZDC[0-7] 
        - VPD[0-7]
        - total bXings


3. Obtain spin patterns

  - Verify we have the fill numbers (fill.txt) and the corresponding
    spin patterns for fills listed in goodruns.dat; if you are unsure: 
    - run "getspinpat" to recreate fill.txt from goodruns.dat
      and download spin patterns from CDEV to ./spinpat
    - append $FMSTXT/fill.txt with ./fill.txt (then cat through uniq to 
      remove duplicated lines.. not necessary, but it's nice to do)
      - cat fill.txt >> $FMSTXT/fill.txt
      - cat $FMSTXT/fill.txt | uniq > $FMSTXT/fill.tmp
      - mv $FMSTXT/fill.{tmp,txt}
    - copy downloaded spinpats to $FMSTXT/spinpat/
  - SEE SPECIAL NOTE BELOW REGARDING F17600
  - run spin_cogging to create spinpat/[fill].spin files from the downloaded CDEV files
    - since STAR spin is opposite source spin and source spin is same as CDEV spin,
      this script implements a sign flip
  - SPECIAL NOTE F17600: I ran a modifier script called mod_17600
    to fix spin pattern for this fill.. see the comments in the
    script for further details; the cdev file is the modified 
    spin pattern for the last 14 runs (*.bad files for the unmodified 
    pattern, where the runs should be omitted from analyses)


4. Accumulate the scaler data into one table: "run accumulate"

  - bunch kicking: you must first generate a list of kicked bunches using "bunch_kicker",
    if one does not exist; the list is in the text file "kicked", with columns
    fill, bx, spinbit
    - use the variable run_randomizer; if it's zero, spinbit equalisation does not run
      (see spinbit equalization section below)
    - bunches which are manually removed are listed in the beginning of the script
    - explanation of algorithm is in the comments

  - execute "accumulate"
    - collects all the datfiles into datfiles/acc.dat, with columns:
      (** and filters out bad part of 17600)
      - run index
      - runnumber
      - fill
      - run time (seconds)
      - bunch crossing (bx)
      - BBC[1-8]
      - ZDC[1-8]
      - VPD[1-4]
      - total bXings
      - blue spin
      - yellow spin

  - accumulate then creates "counts.root," which contains useful trees, using mk_tree.C
  - mk_tree.C reads datfiles/acc.dat file
    - BXING SHIFT CORRECTIONS ARE IMPLEMENTED HERE
    - acc tree: simply the acc.dat table converted into a tree
    - sca tree: restructured tree containing containing branches like 
      - bbc east, bbc west, bbc coincidence 
      - similar entries for zdc and vpd
      - run number, fill number, bunch crossing, spin bit
      - num_runs = number of runs in a fill
      - kicked bunches: 
        certain bunches which are empty according to scalers but filled according to
        cdev are labelled as 'kicked' in the output tree; bXings which are kicked
        are not projected to any distributions in the rellum4 output
  - can run "print_run_table.dat" to print a table containing the columns
    [run index] [run number] [fill index] [fill number]
        

5. Compute the relative luminosity

  - rellum4.C is the analysis script that reads counts.root
  - objects in rdat root file
    - c_spin_pat -- R_spin = same spin / diff spin
    - c_raw_{bbc,zdc,vpd} = raw scaler counts vs var
    - c_acc_{bbc,zdc,vpd} = accidentals corrected scaler counts vs var
    - c_mul_{bbc,zdc,vpd} = multiples corrected scaler counts vs var
    - c_fac_{bbc,zdc,vpd} = correction factor (mult/raw) vs var
    - c_R#_{bbc,zdc,vpd} = relative luminosity vs. var
    - c_mean_R# = mean rellum over EWX (bbc,zdc,vpd on same canvas)
    - c_R#_zdc_minus_vpd = difference between zdc and vpd
    - c_deviation_R#_{bbc,zdc,vpd} = rellum minus mean rellum
    - rate_dep_R#_{bbc*,zdc*,vpd*} = rellum vs. multiples corrected rate
    - c_rate_fac_{bbc*,zdc*,vpd*} = correction factor vs. multiples corrected rate
                                    (tprofile from rate_fac for each spinbit)
  - rellum looping scripts for relative luminosity analysis
    - rellum_all
      - changed often; basically used to run rellum4.C for various independent
        variables etc. 
      - outputs pngs in png_rellum, ready to be copied to protected area to link
        to scalers drupal page
    - rellum_fills
      - runs rellum4.C for all fills separately and output pdfs in subdirectories of
        pdf_bXings_fills; this is for looking at fill dependence of bXing 
        distributions
      - execute "ghost_script" afterward to combine all the pdfs into 
        pdf_bXings_fills/*.pdf
      - this was created to search for the origin of pathologies which cause disagreement
        between R3 and mean R3 (and disagreement between ZDC & VPD?)


6. Combine all the data into a tree to pass to asymmetry analysis

  - run sumTree.C, which builds "sums.root" from "counts.root", which sums the counts
    for each run
    - determines spin pattern types that were collided (see spin pattern
      recognition section below)
    - can now run "nbx_check.C" to test whether the variable "tot_bx" actually
      makes sense with respect to the run time, by plotting tot_bx/(bXing rate) vs.
      run time; the slope of a linear fit to this should equal unity
  - run combineAll.C, which combines sums.root and rdat_i.root into a final tree, which
    can then be passed to asymmetry analysis code


Spinbit Equalizing Running 
--------------------------
 - empty bXings are omitted manually in "bunch_kicker" 
 - number of bXings per spinbit is usually unequal; spinbit equalization randomly 
   removes the minimum number of bXings in order to equalize the number of 
   bXings per spinbit
 - to turn on spinbit equalization, in "bunch_kicker", set run_randomizer=1
   - if run_randomizer!=1, then no bXings other than empty ones will be kicked
 - kicked will now be populated with more bXings to remove; proceed with normal 
   rellum analysis
   - NOTE FOR DRUPAL PAGE: be sure to upload the pngs to protected are in proper directories!
     - scalers2013/png_rellum is not spinbit equalized
     - scalers2013/png_rellum_cleaned_up is the spinbit_equalized (name is old)


Other Useful Scripts
--------------------
 - draw_spin_patterns.C -- draws the spin patterns to "pattern.pdf"
 - draw_fill_vs_run.C -- draws fill index vs. run index
   --> useful for looking for fill structure in plots where 
       run index is the independent variable
 - nbx_check.C plots total bXings / bXing rate vs. run time from sums.root
   --> tau := total bXings / bXing rate (tau should equal run time t)
   --> useful to make sure tot_bx variable makes sense
   --> also looks for runs / fills which have tau != t
 - nbx_check_2.C plots the total # bXings vs. bXing no. for each run 
   into a pdf, called nbx_vs_bxing.pdf --> odd structure?


Making bXing Distributions
--------------------------
- useful for hunting for bad bXings

echo 0 0 0 > kicked
accumulate

rellum_fills with drawLog=0 and zoomIn=0
ghost_script_fills
mv pdf_bXings_fills/*.pdf htmlfiles/pdf_bXings_fills_lin

rellum_fills with drawLog=0 and zoomIn=1
ghost_script_fills
mv pdf_bXings_fills/*.pdf htmlfiles/pdf_bXings_fills_lin_zoom

rellum_fills with drawLog=1 and zoomIn=0
ghost_script_fills
mv pdf_bXings_fills/*.pdf htmlfiles/pdf_bXings_fills_log


Zilong's Scaler Bit Definitions
-------------------------------

- BBC & ZDC -- 3 bits -- [coincidence][west][east]
  - 0 - 000 -- no triggers
  - 1 - 001 -- east
  - 2 - 010 -- west
  - 3 - 011 -- west + east
  - 4 - 100 -- coin
  - 5 - 101 -- coin + east
  - 6 - 110 -- coin + west
  - 7 - 111 -- coin + west + east

- VPD -- 2 bits -- [west][east] (no coincidence)
  - 0 - 00 -- no triggers
  - 1 - 01 -- east
  - 2 - 10 -- west
  - 3 - 11 -- west + east



Spin Pattern Recognition
------------------------

- There are 4 types of spin patterns for Run13:
  - pattern 1: + + - - + + - - 
  - pattern 2: - - + + - - + + 
  - pattern 3: + + - - - - + +
  - pattern 4: - - + + + + - -

- each fill is given an "overall spin pattern no."
  - N := overall spin pattern no.
  - Nb := blue spin pattern no.
  - Ny := yellow spin pattern no.
  - N = 10 * Nb + Ny

- For run13, the following patterns were collided:
  [13, 14, 23, 24, 31, 32, 41, 42]



Deprecated Methods
------------------
  - spinbit_looper
    - loops spinbit script for all fills listed in goodr uns.dat
      - this script will combine the blue and yellow cdev spinpat files
        into one file, named spinpat/[fill#].spin (called *.spin files), with 
        three columns: [bunch] [blue spin] [yellow spin]

