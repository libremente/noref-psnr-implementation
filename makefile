all: 
	gcc no_ref_psnr.c -Wall -lm -o no_ref_psnr
	./no_ref_psnr foreman.xml results 
