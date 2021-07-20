/* Program to do "calculations" on numeric CSV data files.

   Skeleton program written by Alistair Moffat, ammoffat@unimelb.edu.au,
   September 2020, with the intention that it be modified by students
   to add functionality, as required by the assignment specification.

   Student Authorship Declaration:

   (1) I certify that except for the code provided in the initial skeleton
   file, the  program contained in this submission is completely my own
   individual work, except where explicitly noted by further comments that
   provide details otherwise.  I understand that work that has been developed
   by another student, or by me in collaboration with other students, or by
   non-students as a result of request, solicitation, or payment, may not be
   submitted for assessment in this subject.  I understand that submitting for
   assessment work developed by or in collaboration with other students or
   non-students constitutes Academic Misconduct, and may be penalized by mark
   deductions, or by other penalties determined via the University of
   Melbourne Academic Honesty Policy, as described at
   https://academicintegrity.unimelb.edu.au.

   (2) I also certify that I have not provided a copy of this work in either
   softcopy or hardcopy or any other form to any other student, and nor will I
   do so until after the marks are released. I understand that providing my
   work to other students, regardless of my intention or any undertakings made
   to me by that other student, is also Academic Misconduct.

   (3) I further understand that providing a copy of the assignment
   specification to any form of code authoring or assignment tutoring service,
   or drawing the attention of others to such services and code that may have
   been made available via such a service, may be regarded as Student General
   Misconduct (interfering with the teaching activities of the University
   and/or inciting others to commit Academic Misconduct).  I understand that
   an allegation of Student General Misconduct may arise regardless of whether
   or not I personally make use of such solutions or sought benefit from such
   actions.

   Signed by: Harshita Soni. Student ID: 1138784
   Dated:     17 September, 2020

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>

#define MAXCOLS	20	/* maximum number of columns to be handled */
#define MAXROWS	999	/* maximum number of rows to be handled */
#define LABLEN  20	/* maximum length of each column header */
#define LINELEN 100	/* maximum length of command lines */

#define ERROR	(-1)	/* error return value from some functions */

#define O_NOC	'-'	/* the "do nothing" command */
#define O_IND	'i'	/* the "index" command */
#define O_ANA 	'a'	/* the "analyze" command */
#define O_DPY 	'd'	/* the "display" command */
#define O_PLT 	'p'	/* the "plot" command */
#define O_SRT 	's'	/* the "sort" command */

#define CH_COMMA ','	/* comma character */
#define CH_CR    '\r'	/* pesky CR character in DOS-format files */
#define CH_NL    '\n'	/* newline character */

#define TRUE 1
#define FALSE 0
#define EPSILON 1e-6
#define NO_OF_BINS 10
#define MAX_BARS 60
#define HISTOGRAM_UNIT ']'

typedef char head_t[LABLEN+1];   //string type for column headers

typedef double csv_t[MAXROWS][MAXCOLS];   //type for csv data

/****************************************************************/

/* function prototypes */

void get_csv_data(csv_t D, head_t H[],  int *dr,  int *dc, int argc,
                  char *argv[]);
void error_and_exit(char *msg);
void print_prompt(void);
int  get_command(int dc, int *command, int ccols[], int *nccols);
void handle_command(int command, int ccols[], int nccols,
                    csv_t D, head_t H[], int dr, int dc);
void do_index(csv_t D, head_t H[], int dr, int dc, int ccols[], int nccols);

void do_analyze(csv_t D, head_t H[], int dr, int dc,
                int ccols[], int nccols);
void do_display(csv_t D, head_t H[], int dr, int dc,
                int ccols[], int nccols);
void do_sorting(csv_t D, head_t H[], int dr, int dc,
                int ccols[], int nccols);
void do_plotting(csv_t D, head_t H[], int dr, int dc,
                 int ccols[], int nccols);
void get_maxmin_values(csv_t D, int dr, int ccols[], int nccols,
                       double *max, double *min);
void get_limits_and_step(double min, double max, double *lower,
                         double *upper, double *step);
void get_freq_array(csv_t D, int dr, double lower, double upper, double step,
                    int ccols[], int nccols, int freq_array[]);
void print_plot_please(csv_t D, int frq_array[], int n, int ccols[],int nccols,
                       double lower, double upper, double step, int scale);
void print_bar(int frequency, int scale);
void sort_column(csv_t D, int dr, int dc, int c);
void swap_rows(csv_t D, int dr, int dc, int r1, int r2);
void print_display_row(double  row[], int n, int instance);
void copy_row(double from[], double to[], int size);
double get_max(double column[], int ncolumn);
double get_min(double column[], int ncolumn);
double get_avg(double column[], int ncolumn);
double get_median(double column[], int ncolumn);
int is_sorted(double *column, int ncolumn);

/****************************************************************/

/* main program controls all the action */
int
main(int argc, char *argv[]) {

    head_t H[MAXCOLS];	/* labels from the first row in csv file */
    csv_t D;		    /* the csv data stored in a 2d matrix */
    int dr=0, dc=0;		/* number of rows and columns in csv file */
    int ccols[MAXCOLS];
    int nccols;
    int command;

    /* this next is a bit of magic code that you can ignore for
       now, it reads csv data from a file named on the
       commandline and saves it to D, H, dr, and dc
       */
    get_csv_data(D, H, &dr, &dc, argc, argv);

    /* ok, all the input data has been read, safe now to start
       processing commands against it */

    print_prompt();
    while (get_command(dc, &command, ccols, &nccols) != EOF) {
        handle_command(command, ccols, nccols,
                       D, H, dr, dc);
        print_prompt();
    }

    /* all done, so pack up bat and ball and head home */
    printf("\nTa daa!!!\n");
    return 0;
}

/****************************************************************/

/* prints the prompt indicating ready for input */
void
print_prompt(void) {
    printf("> ");
}

/****************************************************************/

/* read a line of input into the array passed as argument
   returns false if there is no input available
   all whitespace characters are removed
   all arguments are checked for validity
   if no arguments, the numbers 0..dc-1 are put into the array
*/
int
get_command(int dc, int *command, int columns[], int *nccols) {
    int i=0, c, col=0;
    char line[LINELEN];
    /* command is in first character position */
    if ((*command=getchar()) == EOF) {
        return EOF;
    }
    /* and now collect the rest of the line, integer by integer,
       sometimes in C you just have to do things the hard way */
    while (((c=getchar())!=EOF) && (c!='\n')) {
        if (isdigit(c)) {
            /* digit contributes to a number */
            line[i++] = c;
        } else if (i!=0)  {
            /* reached end of a number */
            line[i] = '\0';
            columns[col++] = atoi(line);
            /* reset, to collect next number */
            i = 0;
        } else {
            /* just discard it */
        }
    }
    if (i>0) {
        /* reached end of the final number in input line */
        line[i] = '\0';
        columns[col++] = atoi(line);
    }

    if (col==0) {
        /* no column numbers were provided, so generate them */
        for (i=0; i<dc; i++) {
            columns[i] = i;
        }
        *nccols = dc;
        return !EOF;
    }

    /* otherwise, check the one sthat were typed against dc,
       the number of cols in the CSV data that was read */
    for (i=0; i<col; i++) {
        if (columns[i]<0 || columns[i]>=dc) {
            printf("%d is not between 0 and %d\n",
                   columns[i], dc);
            /* and change to "do nothing" command */
            *command = O_NOC;
        }
    }
    /* all good */
    *nccols = col;
    return !EOF;
}

/****************************************************************/

/* this next bit reads the input csv data from a file named on the commandline and saves it into an array of character strings (first line), and into a matrix of doubles (all other lines), using the types defined at the top of the program. 
*/
void
get_csv_data(csv_t D, head_t H[],  int *dr,  int *dc, int argc,
             char *argv[]) {
    FILE *fp;
    int rows=0, cols=0, c, len;
    double num;

    if (argc<2) {
        /* no filename specified */
        error_and_exit("no CSV file named on commandline");
    }
    if (argc>2) {
        /* confusion on command line */
        error_and_exit("too many arguments supplied");
    }
    if ((fp=fopen(argv[1], "r")) == NULL) {
        error_and_exit("cannot open CSV file");
    }

    /* ok, file exists and can be read, next up, first input
       line will be all the headings, need to read them as
       characters and build up the corresponding strings */
    len = 0;
    while ((c=fgetc(fp))!=EOF && (c!=CH_CR) && (c!=CH_NL)) {
        /* process one input character at a time */
        if (c==CH_COMMA) {
            /* previous heading is ended, close it off */
            H[cols][len] = '\0';
            /* and start a new heading */
            cols += 1;
            len = 0;
        } else {
            /* store the character */
            if (len==LABLEN) {
                error_and_exit("a csv heading is too long");
            }
            H[cols][len] = c;
            len++;
        }
    }
    /* and don't forget to close off the last string */
    H[cols][len] = '\0';
    *dc = cols+1;

    /* now to read all of the numbers in, assumption is that the input
       data is properly formatted and error-free, and that every row
       of data has a numeric value provided for every column */
    rows = cols = 0;
    while (fscanf(fp, "%lf", &num) == 1) {
        /* read a number, put it into the matrix */
        if (cols==*dc) {
            /* but first need to start a new row */
            cols = 0;
            rows += 1;
        }
        /* now ok to do the actual assignment... */
        D[rows][cols] = num;
        cols++;
        /* and consume the comma (or newline) that comes straight
           after the number that was just read */
        fgetc(fp);
    }
    /* should be at last column of a row */
    if (cols != *dc) {
        error_and_exit("missing values in input");
    }
    /* and that's it, just a bit of tidying up required now  */
    *dr = rows+1;
    fclose(fp);
    printf("    csv data loaded from %s", argv[1]);
    printf(" (%d rows by %d cols)\n", *dr, *dc);
    return;
}

/****************************************************************/

void
error_and_exit(char *msg) {
    printf("Error: %s\n", msg);
    exit(EXIT_FAILURE);
}

/****************************************************************/

/* the 'i' index command
*/
void
do_index(csv_t D, head_t H[], int dr, int dc,
         int ccols[], int nccols) {
    int i, c;
    printf("\n");
    for (i=0; i<nccols; i++) {
        c = ccols[i];
        printf("    column %2d: %s\n", c, H[c]);
    }
}


/*****************************************************************
******************************************************************

/* this function examines each incoming command and decides what
   to do with it, kind of traffic control, deciding what gets
   called for each command, and which of the arguments it gets
*/
void
handle_command(int command, int ccols[], int nccols,
               csv_t D, head_t H[], int dr, int dc) {
    if (command==O_NOC) {
        /* the null command, just do nothing */
    } else if (command==O_IND) {
        do_index(D, H, dr, dc, ccols, nccols);

    } else if (command==O_ANA) {
        do_analyze(D, H, dr, dc, ccols, nccols);
    } else if (command==O_DPY) {
        do_display(D, H, dr, dc, ccols, nccols);
    } else if (command==O_SRT) {
        do_sorting(D, H, dr, dc, ccols, nccols);
    } else if (command==O_PLT) {
        do_plotting(D, H, dr, dc, ccols, nccols);
        /* and now a last option for things that aren't known */
    } else {
        printf("command '%c' is not recognized"
               " or not implemented yet\n", command);
    }
    return;
}

/* the 'a' index command */
void
do_analyze(csv_t D, head_t H[], int dr, int dc,
           int ccols[], int nccols) {
    double max, min, avg, med;
    /* assign declared arrays to zero to avoid garbage values*/
    double column[MAXROWS]={0};
    int i, j, c, sorted_flag;

    for (i=0; i<nccols; i++) {
        c = ccols[i];
        for (j = 0; j<dr; j++) {
            column[j] = D[j][c];
        }

        max = get_max(column, dr);
        min = get_min(column, dr);
        avg = get_avg(column, dr);
        sorted_flag=FALSE;

        printf("\n");
        printf("         %8s", H[c]);
        /* check if column is sorted and print accordingly */
        if (is_sorted(column, dr)) {
            printf(" (sorted)");
            sorted_flag=TRUE;
        }

        printf("\n");
        printf("    max = %7.1f\n", max);
        printf("    min = %7.1f\n", min);
        printf("    avg = %7.1f\n", avg);
        if (sorted_flag) {
            /* calculate median only if sorted */
            med = get_median(column, dr);
            printf("    med = %7.1f\n", med);
        }
    }
}

/* the 'd' index command */
void
do_display(csv_t D, head_t H[], int dr, int dc,
           int ccols[], int nccols) {
    int i, j, c, instance;
    double first_instance[MAXCOLS]={0}, next_row[MAXCOLS]={0};

    /* printing the column headers according to the specified format */
    printf("\n");
    for (i=nccols-1; i>=0; i--) {
        for (j=i; j>0; j--) {
            printf("        ");
        }
        c = ccols[i];
        printf("%8s\n", H[c]);
    }

    for (i=0; i<nccols; i++) {
        c = ccols[i];
        first_instance[i] = D[0][c];    /* first row */
    }
    instance=1; /* keeps track of the number of times a row is repeated */

    for (i=1; i<dr; i++) {
        /* this loop extracts the next row in the 2d array */
        for (j=0; j<nccols; j++) {
            c = ccols[j];
            next_row[j]=D[i][c];
        }
        /* compare 2 consecutive rows, if not the same, print the first row */
        for (j=0; j<nccols; j++) {
            if (first_instance[j] != next_row[j]) {
                print_display_row(first_instance, nccols, instance);
                /* reset instance. now the next_row array is copied to the
                   first_row and the next entry will be assessed */
                instance=0;
                copy_row(next_row, first_instance, nccols);
            }
        }
        instance++;
    }
    /* don't forget to print the last row! */
    print_display_row(first_instance, nccols, instance);
}

/* the 's' index command */
void do_sorting(csv_t D, head_t H[], int dr, int dc,
                int ccols[], int nccols) {
    int i, c;

    /* sort columns in decreasing order of preference, starting from the last
       entry in ccols. This makes sure that when we reach the primary column
       for sorting, the other columns will already be sorted. E.g, s 3 2 1 is
       dealt with in the following way:
       column 3 is sorted first followed by 2. column 1, which is the primary
       key for sorting is sorted at the last and order is preserved. */
    for (i=nccols-1; i>=0; i--) {
        c = ccols[i];
        sort_column(D, dr, dc, c);
    }
    printf("\n    sorted by: ");
    for (i=0; i<nccols-1; i++) {
        c = ccols[i];
        printf("%s, ", H[c]);
    }
    printf("%s\n", H[ccols[i]]); /* mindful of the commas */
}

/* the 'p' index command */
void do_plotting(csv_t D, head_t H[], int dr, int dc,
                 int ccols[], int nccols) {
    int i, maxfreq;
    double max, min, lower_limit, upper_limit, step;

    /* get the max and min values across all columns to be plotted */
    get_maxmin_values(D, dr, ccols, nccols, &max, &min);
    if (max==min && nccols==1) {
        printf("\nall selected elements are %.1f\n", max);
        return;
    }

    get_limits_and_step(min, max, &lower_limit, &upper_limit, &step);

    /* array of frequencies for each column in each of the 10 bins */
    int freq_array[MAXCOLS*NO_OF_BINS] = {0}, nfreq;

    /* 10(NO_OF_BINS) frequencies for each column */
    nfreq = NO_OF_BINS*nccols;
    get_freq_array(D, dr, lower_limit, upper_limit, step,
                   ccols, nccols, freq_array);

    /* find the maximum value of frequency across all columns */
    maxfreq = freq_array[0];
    for (i = 0; i<nfreq; i++) {
        if (freq_array[i] > maxfreq) {
            maxfreq = freq_array[i];
        }
    }

    /* and then set the scale accordingly. default = 1 */
    int scale=1; double ratio=0.0;
    if (maxfreq > MAX_BARS) {
        ratio = ((double)maxfreq)/MAX_BARS;
        scale = ceil(ratio);
    }
    print_plot_please(D, freq_array, nfreq, ccols, nccols,
                      lower_limit, upper_limit, step, scale);
}


/* Finds the maximum and minimum value across all columns passed in ccols[]
   Helper function to the do_plotting() function for 'p' command */
void get_maxmin_values(csv_t D, int dr, int ccols[], int nccols,
                       double *max, double *min) {
    int i, j, c;
    double localmax, localmin, column[MAXROWS]={0};
    c = ccols[0];
    for (j=0; j<dr; j++) {
        column[j] = D[j][c];
    }
    *max = get_max(column, dr);
    *min = get_min(column, dr);

    for (i=1; i<nccols; i++) {
        c = ccols[i];
        for (j=0; j<dr; j++) {
            column[j] = D[j][c];
        }
        localmax = get_max(column, dr);  /* max of the particular column */
        localmin = get_min(column, dr);
        if (localmax > *max)
            *max = localmax;
        if (localmin < *min)
            *min = localmin;
    }
}

/* Calculates the step, upper and lower limits for plotting in 'p' command*/
void get_limits_and_step(double min, double max, double *lower,
                         double *upper, double *step) {
    *lower = min - EPSILON;
    *upper = max + EPSILON;
    *step = (*upper - *lower) / NO_OF_BINS;
}

/* Calculates the frequency of each column in each bin in ccols[].
   Size of freq_array[] will be (no. of bins)*(no. of columns).
   Helper function to the do_plotting() function for 'p' command */
void get_freq_array(csv_t D, int dr, double lower, double upper, double step,
                    int ccols[], int nccols, int freq_array[]) {
    int i, j, c, frequency, n=0;
    double curr, next;
    printf("\n");
    curr=lower;
    while(curr<upper) {
        for (i=0; i<nccols; i++) {
            next = curr+step; /* upper bound for current bin */
            c = ccols[i];
            frequency = 0;
            for (j = 0; j<dr; j++) {
                if (curr <= D[j][c] && D[j][c]<next) {
                    frequency+=1;
                }
            }
            freq_array[n++] = frequency;
        }
        curr+=step;
    }
}

/* helper function to print the plot once we have the bins and frequencies */
void print_plot_please(csv_t D, int frq_array[], int n, int ccols[],int nccols,
                       double lower, double upper, double step, int scale) {
    int i, j, frequency;
    double curr;
    curr=lower; /* lower end of the starting bin */

    /* iterate through the freq_array. frequency values of a column for each
       bin will be located 'nccols' index apart */
    for (j=0; j<n; j+=nccols) {
        printf("    %7.1f +\n", curr);
        for (i=0; i<nccols; i++){  /* print and plot each column in ccols[] */
            printf("         %2d |", ccols[i]);
            frequency = frq_array[i+j];
            if (frequency) {
                print_bar(frequency, scale);
            }
            printf("\n");
        }
        curr+=step;  /* move onto next bin */
    }
    printf("    %7.1f +\n", upper);  /* upper end of the last bin */
    printf("    scale = %d\n", scale);
}

/* Prints a bar of the histogram. Helper func for the 'p' command */
void print_bar(int frequency, int scale) {
    int i, no_of_bars = ceil((double)frequency/(double)scale);
    for(i=0; i<no_of_bars; i++) {
        printf("%c", HISTOGRAM_UNIT);
    }
}

/* Sorts the 2d arrray in ascending order on the basis of a column.
   Column number = int c */
void sort_column(csv_t D, int dr, int dc, int c) {
    int i, j;
    /* modified version of insertion sort code taken from:
       https://people.eng.unimelb.edu.au/ammoffat/ppsaa/c/insertionsort.c */
    for (i = 1; i<dr; i++) {
        for (j = i-1; j>=0; j--) {
            if (D[j+1][c] < D[j][c]) {
                /* swap entire rows */
                swap_rows(D, dr, dc, j+1, j);
            }
        }
    }
}

/* Swaps two given rows using memcpy() function from <string.h> */
void swap_rows(csv_t D, int dr, int dc, int r1, int r2) {
    double temp[MAXCOLS];
    memcpy(temp, D[r1], dc*sizeof(double ));
    memcpy(D[r1], D[r2], dc*sizeof(double ));
    memcpy(D[r2], temp, dc*sizeof(double ));
}

/* Prints an output row as per specified format for 'd' command.
   Helper function to the do_display() function */
void print_display_row(double  row[], int n, int instance) {
    int i;
    for (i=0; i<n; i++) {
        printf(" %7.1f", row[i]);
    }
    if (instance == 1) {
        printf("    (%2d instance)\n", instance);
    } else {
        printf("    (%2d instances)\n", instance);
    }
}

/* Copies an array(from) to another(to) using memcpy() from <string.h> */
void copy_row(double from[], double to[], int size) {
    memcpy(to, from, size*sizeof(double));
}

/* Finds the maximum value of the array passed (column) */
double get_max(double column[], int ncolumn) {
    double max = column[0];
    int i;
    for (i=1; i<ncolumn; i++) {
        if (column[i] > max) {
            max = column[i];
        }
    }
    return max;
}

/* Finds the minimum value of the array passed (column) */
double get_min(double column[], int ncolumn) {
    double min = column[0];
    int i;
    for (i=1; i<ncolumn; i++) {
        if (column[i] < min) {
            min = column[i];
        }
    }
    return min;
}

/* Calculates the average of the array passed (column) */
double get_avg(double column[], int ncolumn) {
    double avg, sum=0.0;
    int i;
    for (i=0; i<ncolumn; i++) {
        sum += column[i];
    }
    avg = sum/ncolumn;
    return avg;
}

/* Calculates the median of the array passed (column) */
double get_median(double column[], int ncolumn) {
    double median;
    if (ncolumn % 2 == 1) {
        /* need to do -1 rather than +1 because arrays are zero-offset */
        median = (column[(ncolumn/2)-1] + column[ncolumn/2])/2;
    } else {
        median = column[(ncolumn - 1)/2];
    }
    return median;
}

/* Recursively checks whether the array(column) is sorted in ascending order */
int is_sorted(double *column, int ncolumn) {
    if (ncolumn<=1) {
        return TRUE;
    }
    if (column[ncolumn-1] < column[ncolumn-2]) {
        return FALSE;
    }
    return is_sorted(column, (ncolumn-1));
}

/* ALGORITHMS ARE FUN! */
