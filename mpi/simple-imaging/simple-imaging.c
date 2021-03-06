
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include <stdint.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <complex.h>
#include <fftw3.h>

#include "grid.h"

#include <mpi.h>

/* ****************************************************************************************************************** */
/* Parallel MPI version for simple imaging, this version assumes that the visibilities have already been distributed  */
/* by frequency and time and each node reads a file which cointains the a subset of the visibilities                  */
/* The visibility file is the same file for all nodes at this point 						      */
/* TODO: Generate a realistic input with Oscar so that the image makes sense					      */
/* ****************************************************************************************************************** */


int main(int argc, char *argv[]) {

    // Read parameters
    static struct option options[] =
      {
        {"theta",   required_argument, 0, 't' },
        {"lambda",  required_argument, 0, 'l' },
        {"wkern",   optional_argument, 0, 'w' },
        {"akern",   optional_argument, 0, 'a' },
        {"grid",    optional_argument, 0, 'g' },
        {"image",   optional_argument, 0, 'i' },
        {"min-bl",  optional_argument, 0, 'b' },
        {"max-bl",  optional_argument, 0, 'B' },
        {0, 0, 0, 0}
      };
    int option_index = 0;
    //double theta = 0, lambda = 0;
	double theta = 0.002; //0.08; // Width of the Field of view as directional cosines (aproximately radians)  (0.002 in the notebook aaf for generated vis)
	double lambda = 15000; // 300000;  //  Width of the uv-plane (in wavelengths). Controls resolution of the images. (15000 in the notebook aaf)
    char *wkern_file = NULL, *akern_file = NULL,
         *grid_file = NULL, *image_file = NULL;
    double bl_min = DBL_MIN, bl_max = DBL_MAX;
    int c; int invalid = 0;
    while ((c = getopt_long(argc, argv, ":", options, &option_index)) != -1) {
        switch(c) {
        case 't': theta = atof(optarg); break;
        case 'l': lambda = atof(optarg); break;
        case 'w': wkern_file = optarg; break;
        case 'a': akern_file = optarg; break;
        case 'g': grid_file = optarg; break;
        case 'i': image_file = optarg; break;
        case 'b': bl_min = atof(optarg); break;
        case 'B': bl_max = atof(optarg); break;
        default: invalid = 1; break;
        }
    }

	//MPI
	int myid, numprocs;
  	MPI_Status status;
	int provided;

  	//MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  	//assert(MPI_THREAD_MULTIPLE == provided);
  	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  	MPI_Comm_rank(MPI_COMM_WORLD, &myid);



    // Check grid parameters
    int grid_size = (int)(theta * lambda);
    size_t grid_byte_size = grid_size * grid_size * sizeof(double complex);
    if (grid_size <= 0) {
        fprintf(stderr, "Invalid grid configuration!\n");
        invalid = 1;
    }

#ifdef READ_VIS
    // Must have an input file
    const char *vis_file = 0;
    if (optind + 1 == argc) {
        vis_file = argv[optind];
    } else {
        printf("Please supply a visibility input file!\n");
        invalid = 1;
    }
#endif

    if (invalid) {
        printf("usage: %s --theta=THETA --lambda=LAM [--grid=GRID]\n", argv[0]);
        printf("              [--image=IMAGE] [--wkern=WKERN] [--akern=AKERN]\n");
        printf("              [--min-bl=MIN_BL] [--max-bl=MAX_BL]\n");
        printf("              INPUT\n");
        printf("\n");
        printf("optional arguments:\n");
        printf("  --theta=THETA         Field of view size (in radians)\n");
        printf("  --lambda=LAM          uv grid size (in wavelenghts)\n");
        printf("  --grid=GRID           grid output file\n");
        printf("  --image=IMAGE         image output file\n");
        printf("  --wkern=WKERN         w-kernel file to use for w-projection\n");
        printf("  --akern=AKERN         A-kernel file to use for w-projection\n");
        printf("  --min-bl=MIN_BL       Minimum baseline length to consider (in km)\n");
        printf("  --max-bl=MAX_BL       Maximum baseline length to consider (in km)\n");
        printf("positional arguments:\n");
        printf("  input                 input visibilities\n");
        return 1;
    }

    // Intialise HDF5
    init_dtype_cpx();

    // Open files
    struct vis_data vis;
    struct w_kernel_data wkern;
    struct a_kernel_data akern;
    int grid_fd = -1, image_fd = -1;

	// FIXME: All procs read same file
	// Open files
	// Read visibilities from file
#ifdef READ_VIS
    if (load_vis(vis_file, &vis, bl_min, bl_max)) {
        return 1;
    }
#else
	// Generate visibilities
	if(generate_vis(NULL, &vis, bl_min, bl_max)){
		printf("Error generating visibilities\n");
		exit(0);
	}	

#endif
    if (wkern_file) {
        if (load_wkern(wkern_file, theta, &wkern)) {
            return 1;
        }
    }
    if (wkern_file && akern_file) {
        if (load_akern(akern_file, theta, &akern)) {
            return 1;
        }
    }
    if (grid_file) {
        grid_fd = open(grid_file, O_CREAT | O_TRUNC | O_WRONLY, 0777);
        if (grid_fd == -1) {
            perror("Failed to open grid file");
            return 1;
        }
    }
    if (image_file) {
        image_fd = open(image_file, O_CREAT | O_TRUNC | O_WRONLY, 0777);
        if (image_fd == -1) {
            perror("Failed to open image file");
            return 1;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);  // For the timing
    double t,t1,t2;
    if(myid ==0){ 
	t1 = wtime();
    }

    // Both master and worker work on a subset of visibilities 

    // Allocate grid
    printf("\nGrid size:    %d x %d (%.2f GB)\n", grid_size, grid_size, (double)(grid_byte_size)/1000000000);
    double complex *uvgrid = (double complex *)calloc(grid_byte_size, 1);

    // Simple uniform weight (we re-use the grid to save an allocation)
    printf("Weighting...\n");
    weight((unsigned int *)uvgrid, grid_size, theta, &vis);
    memset(uvgrid, 0, grid_size * grid_size * sizeof(unsigned int));

    // Set up performance counters
//    struct perf_counters counters;
//    open_perf_counters(&counters);



    uint64_t flops = 0, mem = 0;
    if (!wkern_file) {
        printf("Gridder: Simple imaging\n");
  //      enable_perf_counters(&counters);
        flops = grid_simple(uvgrid, grid_size, theta, &vis);
  //      disable_perf_counters(&counters);
        // Assuming 0.5 flop/B
        mem = flops * 2;
    } else if (!akern_file) {
        printf("Gridder: W-projection\n");
  //      enable_perf_counters(&counters);
        flops = grid_wprojection(uvgrid, grid_size, theta, &vis, &wkern);
  //      disable_perf_counters(&counters);
        // Assuming 10 flop/B
        mem = flops / 10;
    } else {
        printf("Gridder: AW-projection (ignoring convolution)\n");

        // We need to chunk the convolution so we don't run out of
        // memory...
        const int bl_chunk = 1000;
        int bl_min;
        for (bl_min = 0; bl_min < vis.bl_count; bl_min+=bl_chunk) {
            int bl_max = bl_min + bl_chunk;
            if (bl_max > vis.bl_count) bl_max = vis.bl_count;

            printf("\"Convolving\" %d-%d...\n", bl_min, bl_max-1);
            int bl;
            for (bl = bl_min; bl < bl_max; bl++) {
                convolve_aw_kernels(&vis.bl[bl], &wkern, &akern);
            }

            // Do convolution
    //        enable_perf_counters(&counters);
            flops += grid_awprojection(uvgrid, grid_size, theta, &vis, &wkern, &akern, bl_min, bl_max);
    //        disable_perf_counters(&counters);

            // Free convolved kernels
            for (bl = bl_min; bl < bl_max; bl++) {
                free(vis.bl[bl].awkern);
                vis.bl[bl].awkern = NULL;
            }
        }

        // Assuming 0.5 flop/B
        mem = flops * 2;
    }

    // Show performance counters after gridding
 //   printf("\nCounters:\n");
 //   print_perf_counters(&counters, flops, mem);

    // Make hermitian
    printf("\nMake hermitian...\n");
    make_hermitian(uvgrid, grid_size);

    // Each process Writes his grid
    if (grid_fd != -1) {
        printf("Write grid...\n");
        int i;
        for (i = 0; i < grid_size; i++) {
            write(grid_fd, uvgrid + i * grid_size, grid_byte_size / grid_size);
        }
        close(grid_fd);
    }
        printf("FFT...\n");

        // First shift zero frequency
        fft_shift(uvgrid, grid_size);

        // Do DFT. Complex-to-complex to keep with numpy (TODO: optimize)
        fftw_plan plan;
        plan = fftw_plan_dft_2d(grid_size, grid_size, uvgrid, uvgrid, -1, FFTW_ESTIMATE);
        fftw_execute_dft(plan, uvgrid, uvgrid);

        // Shift zero frequency back into centre
        fft_shift(uvgrid, grid_size);

    double complex *imagegrid=NULL;    
    if(myid == 0){
	imagegrid= (double complex *)calloc(grid_byte_size, 1);
    }

    MPI_Reduce(uvgrid,imagegrid,grid_size*grid_size,MPI_C_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD);
    // Get the image by getting the real part of the uvgrid
    // and MPI_Reduction to accumulate images
    // or accumulate uvgrids and get the image
    // FIXME We get all the grid, see if we can reduce just the image!!!

    // MPI_Reduction to accumulate all images

    if(myid ==0){  // Master prints image
      if (image_fd != -1) {
        // Write real part to disk
        printf("Write image...\n");
        int i;
        double *row = malloc(sizeof(double) * grid_size);
        for (i = 0; i < grid_size; i++) {
            int j;
            for (j = 0; j < grid_size; j++) {
                row[j] = creal(uvgrid[i*grid_size+j]);
            }
            write(image_fd, row, sizeof(double) * grid_size);
        }
        close(image_fd);
     }
     t2=wtime();
     fprintf(stdout, "Time: %04.3f sec", t2-t1);
     fprintf(stdout, "(%3.3f GFlop => %6.2f MFlop/s)\n",
            flops/1000000000.0,
            flops/(t2-t1)/1000000);

    }
    MPI_Finalize();
    return 0;
}
