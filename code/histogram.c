/* Functions to calculate the histogram of PWM scores on a set of sequences
 * */

#include <float.h>
#include "options.h"
#include "alloc-error.h"

static double binwidth;                             /* Width of a bin */
static double num_lmers = 0.0;                      /* counts the number of
                                                       sequences 
                                                       encountered */
static double *histogram;                           /* A double array to 
                                                       hold frequency 
                                                       histogram */

static double *bin_starts;                          /* Lower bounds of bin
                                                       ranges */

static double max_double_int;                       /* the maximum integer
                                                       such that it, and
                                                       everything less than
                                                       is can be
                                                       represented by a
                                                       double           */


void init_histogram();                         /* Initialize 
                                                       histogram */
void update_histogram(double score);                /* Count the score */
void write_histogram();                  /* Write out the
                                                       histogram to file
                                                       fname */

void free_histogram();                              /* Free the histogram
                                                                        */

void init_uniform_bins();                           /* Divide the score
                                                       range into 
                                                       equally-sized bins
                                                                        */
int find_bin(double s);                             /* Find a score's bin
                                                                        */
                                                                        

/* set the number of bins and allocate memory to the histogram */
void init_histogram()
{

    int j;

    /* if hist_numbins is zero, divide the score range into 0.1-wide bins */
    if (hist_numbins == 0)
        hist_numbins = (long int) 10*(ceil(Max_exact_score) - 
                                        floor(Min_exact_score));

    /* Allocate histogram array */
    histogram = (double *) calloc_error(hist_numbins,
                                        sizeof(double),
                                        "histogram",
                                        "init_histogram()");

    /* Allocate bin_starts array */
    bin_starts = (double *) calloc_error(hist_numbins,
                                        sizeof(double),
                                        "bin_starts",
                                        "init_histogram()");

    for (j = 0; j < hist_numbins; j++) {
        histogram[j] = 0.0;
        bin_starts[j] = 0.0;
    }    

    /* populate bin_starts by calling the appropriate method */
    if (!strcmp(binmethod, "uniform"))
        init_uniform_bins();
    else {
        fprintf(stderr, "init_histogram(): Unknown binning method " \
                                                    "%s\n", binmethod);
        exit(1);
    }    

    /* calculate the maximum integer representable by a double less one*/
    max_double_int = pow(FLT_RADIX, DBL_MANT_DIG) - 1.0;

}                                        

/* Free the histogram */
void free_histogram()
{

    free_error(histogram, "histogram", "free_histogram");
    free_error(bin_starts, "bin_starts", "free_histogram");

}

/* Divide the score range into equally sized bins and populate bin_starts
 * array. By choosing integer start and end, this gives bins of width 0.1
 * when number of bins has not been specified */
void init_uniform_bins()
{

   int j;
   double width;

   width = (ceil(Max_exact_score) - floor(Min_exact_score))/hist_numbins;
   for (j = 0; j < hist_numbins; j++) 
       bin_starts[j] = floor(Min_exact_score) + width*j;  

}

/* count the occurance of score and update number of sequences */
void update_histogram(double score)
{

    int binindex;

    binindex = find_bin(score);

    /* DEBUGGING code below, should be handled by a DEBUG flag in future */
    /*
    if ((binindex == 0) && (score >= bin_starts[1]))
        printf("Start [%f, %f, %f]\n", bin_starts[0], score, bin_starts[1]);
    else if ((binindex == hist_numbins - 1) &&
                        (score < bin_starts[hist_numbins - 1]))
        printf("End [%f, %f, Inf]\n", bin_starts[hist_numbins - 1], score);
    else if ( ( (binindex != hist_numbins - 1) && (binindex != 0) ) && 
                        ( (score < bin_starts[binindex]) || 
                            (score >= bin_starts[binindex + 1]) ) )
        printf("Middle [%f, %f, %f]\n", bin_starts[binindex], score, bin_starts[binindex + 1]);
    */    

    histogram[binindex] += 1.0;
    num_lmers += 1.0;

    if (num_lmers > max_double_int) {

        fprintf(stderr, "update_histogram(): More L-mers than " \
                                    "can be represented (%.0f)!\n",
                                                    max_double_int);
        exit(1);

    }    
    
}

/* Do a binary search and return the index of the bin for a given score */
int find_bin(double s)
{

    long int n, midpoint, left, right;

    left = 0;
    right = hist_numbins - 1;
    n = hist_numbins;

    while (n > 1) {

        midpoint = (n % 2 == 0) ? left + n/2 : left + n/2 + 1;

        if (s >= bin_starts[(int) midpoint]) 
            left = midpoint;
        else
            right = midpoint - 1;

        n = right - left + 1;
    }

    return (int) left;

}
    
/* If a filename is specified, write to it, otherwise, write to stdout */
void write_histogram()    
{

    FILE *fp;
    int j;

    if (Hist_file == (char *) NULL)
        fp = stdout;
    else
        fp = fopen(Hist_file, "w");

    /* Write header */
    fprintf(fp, "#Frequency histogram for PWM %s running on " \
                "sequence file %s\n", Mat_file, Seq_file);

    /* Write out the frequency histogram */
    fprintf(fp, "#Bin_start\t\tFrequency\n");
    for (j = 0; j < hist_numbins; j++)
        fprintf(fp, "%.6f\t\t%.0f\n", bin_starts[j], 
                                      histogram[j]);

    if (Hist_file == (char *) NULL)
        fclose(fp);

}
