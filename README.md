# Nesting_function <br>
  2020.12 MANO Poobodin <br>
=============================================== <br>
 Fermi surface calculating program <br>
 & Nesting function calculating program
 Read eigen values from Quantum Espresso / VASP output
 READED files [ Quantum Espresso ]
 1. prefix.band       < from bands.x   >
 2. prefix.band.proj  < from projwfc.x >
 READED files [ VASP ]
 1. PROCAR
===============================================

 &input_plot <br>
 file_read     : prefix of files to be read <br>
 file_k        : name of Fermi surface output file <br>
 mode          : 'qe' , 'vasp' <br>
 ldecom        : 0 for non-decommposed,
                 1 for decomposed Fermi surface
 lsuscep       : 0 for no lindhard or nesting,
                 1 for including lindhard or nesting
 ef            : Fermi energy
 degauss       : degauss for constaining calculation to Fermi level
 thresh        : exclude small value from output-weight
 n_con         : number of orbital to be considered in decomposed bands
 /
 LIST OF DECOMPOSED BANDS
 &suscep_q
 mode          : 'lindhard' or 'nesting'
 k_resoved     : 0 : NO , 1 : YES
 ibnd_start    : start band in band loop
 ibnd_final    : end band in band loop (This makes calculation faster)
 ndegauss      : number of degauss used in calculation
 degauss_init  : initial value for degauss
 delta_gauss   : interval for each degauss
 window        : energy window at Fermi level
 divide_k      : use subset of total k-point
 divide_q      : use subset of total q-point. Useful in convergence test
                 and in speeding up to calculation
 file_q        : output file for lindhard/nesting function
 file_conv     : output file for convergence test
 omega         : coefficient for imaginary part
 /
