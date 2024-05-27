The parameter file "f0.05meV_DnC_repeat_b0.2.param" can be used to generate the quantum dot spetra under strong driving depicted in Fig. 4 (timings in Fig. 5) in the article [Phys. Rev. X 14, 011010 (2024)].
To this end, run

ACE f0.05meV_DnC_repeat_b0.2.param, remove data points with negative time, and Fourier transform the data, e.g., using the tool

FT_outfile -infile f0.05meV_DnC_repeat_b0.2.out -outfile f0.05meV_DnC_repeat_b0.2.fft -sign -1 -t_start 0 -subtract_final true



The files "compare_DnC_b0.2.param", "compare_DnC_b1.param", and "compare_sequential.param", can be used to compare computation times. These correspond to the slice at n=1024 time steps in Fig. 5. Of course, the computation times strongly depend on the CPU used.



