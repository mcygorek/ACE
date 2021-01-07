This set of examples describes the resonant-level model in the ACE article:

First, run the follwing example:

> ACE 01_resonantlevelmodel_N2_degenerate.param 

This describes a single site coupled to two other sites at the same energy 
(Fig 2a in the ACE article)

Column 1 and 2 of the output file '01_resonantlevelmodel_N2_degenerate.out'
will contain the time and the site occupations.



The curves for Fig 2b in the ACE article can be generated similarly by running:

> ACE 02_resonantlevelmodel_N4.param 

> ACE 03_resonantlevelmodel_N10.param 

> ACE 04_resonantlevelmodel_N100.param

Note that the last one requires some calculation time and memory.



