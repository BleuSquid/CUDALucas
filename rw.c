/*!!time & date stamp error messages? */
/*!!Other input formats:
  M( exponent )letter: k=value [k2=value]
  M( exponent )Z [no more input for exponent]
  M( M( exponent ) )[...]
  all factors 1
  all factors 0
  [other lines that change flags given on command line?]
 */

static const char Options[] =
"-\t\tread stdin (default if no files or exponents given)\n"
"-b\t\tuse binary Lucas-Lehmer checkpoint files\n"
"-c#\t\tcheckpoint LL test every # iterations\n"
"-d\t\tduplicate output to stderr\n"
"-d-\t\tduplicate output to stderr\n"
"-dstderr\tduplicate output to stderr\n"
"-dstdout\tduplicate output to stdout\n"
"-d<fn>\t\tappend duplicate output to file <fn>\n"
"-h\t\tuse human-readable LL checkpoint files (fftlucas only)\n"
"-o\t\tsend output to stderr\n"
"-o-\t\tsend output to stdout\n"
"-ostderr\tsend output to stderr\n"
"-ostdout\tsend output to stdout\n"
"-o<fn>\t\tappend output to file <fn>\n"
"-r# #\t\tskip exponents outside range from first # to second #\n"
"-s\t\tsmall test (mostly for the Makefile to test factorers)\n"
"-t\t\tprint CPU time, wall clock time, etc. if available\n"
"-v\t\tprint brief version information\n"
"<fn>\t\tread from file <fn> (Mersenne exponents, 'extract -2'\n"
"\t\toutput, http://www.garlic.com/~wedgingt/mersfmt.html formats)\n"
"<p>\t\ttest Mersenne exponent <p>\n"
"\t\t(<p> will be tried as a file first)\n";

#include "setup.h"
#include "balance.h"
#include "rw.h"

#ifdef _MSC_VER
#include "timeval.h"
#else
#include <sys/time.h>
#endif

#ifdef _MSC_VER
#define UNLINK _unlink
#else
#define UNLINK unlink
#endif

static timeval start_time;

#if defined(USE_RUSAGE)
# include <sys/resource.h>
int getrusage (int who, struct rusage *usage);
#endif

/* Checkpoints are the save files for the Lucas-Lehmer test programs, */
/*  mersenne1, mersenne2, and fftlucas */

UL chkpnt_iterations = 0;

int device_number = 0; //msft

/* checkpoint file formats: human readable and machine-specific binary with magic number ala file(1) */
#define CHKPNT_ASCII 0
#define CHKPNT_BINARY 1
#define MAGIC_NUMBER (0x76757473)

static int chkpnt_fmt = CHKPNT_BINARY;
static char chkpnt_cfn[CHKPNT_FILENAME_LENGTH];
static char chkpnt_tfn[CHKPNT_FILENAME_LENGTH];
static EXPONENT start_range = (UL)31, /* smallest exponent that can be tested */
  stop_range = (~(UL)0)/2; /* this is probably the largest */


static int print_times = 0; /* print wall clock and CPU times to stderr? -t flag */

void print_time_from_seconds(int sec)
{
  if (sec > 3600)
   {
    printf("%d", sec/3600);
    sec %= 3600;
    printf(":%02d", sec/60);
   }
  else
    printf("%d", sec/60);
  sec %= 60;
  printf(":%02d", sec);
}

void print_time (int iterations, int current, int total)
{
  timeval end_time;
  unsigned long diff, diff1;
#if defined(RUSAGE_SELF)
  struct rusage ru;
#endif

  if (!print_times) return;

#if defined(RUSAGE_SELF)
  getrusage(RUSAGE_SELF, &ru);
  printf("%ld.%02ld user %ld.%02ld sys ", (long)ru.ru_utime.tv_sec, (long)ru.ru_utime.tv_usec/10000, (long)ru.ru_stime.tv_sec, (long)ru.ru_stime.tv_usec/10000);
#endif
  gettimeofday(&end_time, NULL);

  diff = end_time.tv_sec - start_time.tv_sec;
  diff1 = 1000000*diff + end_time.tv_usec - start_time.tv_usec;
  start_time = end_time;
  print_time_from_seconds(diff);
  printf(" real");
  if(iterations)
    printf(", %3.4f ms/iter", diff1/1000.0/iterations);
  if(total)
  {
    diff = (long)((total - current)/iterations*(diff1/1e6));
    printf(", ETA ");
    print_time_from_seconds(diff);
  }
  return;
}

/*!!this is a function to make it easier to modify for differing operating systems */
const char *archive_name()
{
  return("mersarch.txt");
}

FILE *open_archive()
{
  return(fopen(archive_name(), "a"));
}

int read_check_point(FILE * infp, UL * q, UL * n, UL * j, double * err, double **x)
{
  /* check for binary checkpoint format */
  UL tmp2;
  if(fread(&tmp2, sizeof(tmp2), 1, infp) == 1
     && tmp2 == MAGIC_NUMBER
     && fread(q, sizeof(*q), 1, infp) == 1
     && fread(n, sizeof(*n), 1, infp) == 1
     && *n < ((1L << 31) + 1L)
     && fread(j, sizeof(*j), 1, infp) == 1
     && fread(err, sizeof(*err), 1, infp) == 1
     && ((*x) == NULL || (free((char *)(*x)), 1), (((*x) = (double *)calloc((size_t)(*n + *n), sizeof(double))) != NULL)
     && fread(*x, sizeof(x[0][0]), (size_t)*n, infp) == *n)
     && (*n) > ((double)(*q))/ROUNDED_HALF_MANTISSA + 3
     && *q > (UL)30
     && *q % 2 == 1
     && *j < *q
     && *j > 0
     && (*err < 0.3 || (*q > (UL)10000 && *err < 0.49))
  )
  {
    for (tmp2 = *n; (tmp2 & 1) == 0; tmp2 >>= 1);

    sprintf(chkpnt_cfn, "c" PRINTF_FMT_UL, *q);
    sprintf(chkpnt_tfn, "t" PRINTF_FMT_UL, *q);

    if (tmp2 == 1 || tmp2 == 3 || tmp2 == 5 || tmp2 == 7)
    {
      fclose(infp);
      return 2;
    }

    fprintf(stderr, "%s: inconsistent binary checkpoint; bad FFT run length of " PRINTF_FMT_UL "\n", program_name, *n);
    return 1;
  }

  if (tmp2 == MAGIC_NUMBER || ferror(infp) || feof(infp))
  {
    fprintf(stderr, "%s: inconsistent or incomplete binary checkpoint\n", program_name);
    return 1;
  }
  return 0;
}

/* Get some input and/or parse the command line (since exponents can be on the command line) */
/*  return values: */
/*  0: no more input; 1: trouble; 2: resuming from a partial (RI) result; 3: other data for exponent */
int input( int argc, char **argv
         , EXPONENT *q, UL *n, UL *j, double *err, double **x, EXPONENT last, char *M
         , FILE **infp, FILE **outfp, FILE **dupfp, const char *version_info)
{
  if (q == NULL || M == NULL || infp == NULL || outfp == NULL || dupfp == NULL) return 1;

  static bool first_time_here = true; /* first call to input() in this program? */
  char *tmpstr = NULL;
  if (first_time_here)
  {
    gettimeofday(&start_time, NULL);
    if (((tmpstr = getenv("FFTLUCAS_CHECKPOINT")) != NULL) && (sscanf(tmpstr, SCANF_FMT_UL, &chkpnt_iterations) != 1 || chkpnt_iterations < (UL)5000))
    {
      if (!isdigit(tmpstr[0]))
        fprintf(stderr, "%s: non-numeric FFTLUCAS_CHECKPOINT; using 5000\n", program_name);
      chkpnt_iterations = (UL)5000;
    }
  }  

  int sBL, sD; /* sizeof(BIG_LONG), sizeof(double), and */
  long RHM;     /* ROUNDED_HALF_MANTISSA as read from input */
  UL tmp2;
  FILE *tmpfp = NULL;
  int partial; /* are we resuming a partially computed LL test? */
  int i = EOF;
  /* these last three are needed by LL tests only to skip parts of the input appropriately */
  char str[1]; /* for reading the ending 's' of 'trials' from input (usually the server) */

get_next_q: /* all but one goto to here are just before a return if (*q < start_range || *q > stop_range); */
            /*  the single exception is reading a DATABASE file that has exponents that do not fit in a UL */

  UL tmp = 0;
  static int arg = 0; /* subscript of argv; will be incremented immediately */
  if ((*infp == NULL || feof(*infp)) && argv != NULL && argv[arg] != NULL)
    arg++;

  while (argv != NULL && argv[arg] != NULL && argv[arg][0] == '\0')
    arg++;

  if ((argv == NULL || argv[arg] == NULL) && (*infp == NULL || feof(*infp)))
  {
    if (first_time_here)
    {
      *infp  = stdin;
      *outfp = stdout;
    }
    else
      return 0; /* nothing else to read; must be done; only return(0) when no -n is given */
  }

  while (arg < argc && argv != NULL && argv[arg] != NULL && ((*infp == NULL && tmp == (UL)0) || (*infp != NULL && feof(*infp)) || first_time_here))
  {
    if (*infp != NULL && *infp != stdin && feof(*infp))
    {
      fclose(*infp);
      *infp = NULL;
    }
    if (*outfp != NULL && *outfp != stdout && *outfp != stderr && ferror(*outfp))
    {
      fclose(*outfp);
      *outfp = NULL;
    }
    switch (argv[arg][0])
    {
      case '\0': /* empty argument: pretend it's a lone dash */
        *infp = stdin;
        first_time_here = false;
        if (arg > 1)
          fprintf(stderr, "Warning: reading from stdin due to NULL argument or no file arguments\n");
        break;
      case '-': /* leading dash: some sort of option */
        /* set tmpstr to point at the option's likely argument to simplify things later */
        if (argv[arg][1] != '\0') /* could skip this except: */
        {
          if (argv[arg][2] != '\0') /* can't check this unless prior condition is true */
            tmpstr = argv[arg] + 2;
          else if (arg + 1 < argc)
            tmpstr = argv[arg + 1];
          else
            tmpstr = NULL;
        }
        switch (argv[arg][1])
        {
          case '\0': /* just a dash: read stdin */
            *infp = stdin;
            first_time_here = 0;
            /*!!kludge to undo increment just after this switch so increment at EOF will be correct */
            /*!! and to match argv[arg] being the name of the open file in the named file case */
            arg--;
            break;
          case 'b': /* binary checkpoint format */
            chkpnt_fmt = CHKPNT_BINARY;
            break;
          case 'c': /* check point iterations */
            if (tmpstr == NULL || sscanf(tmpstr, SCANF_FMT_UL, &chkpnt_iterations) != 1)
              chkpnt_iterations = (UL)5000;
            else
              if (tmpstr != NULL && tmpstr == argv[arg + 1])
                arg++;
            break;
           case 'D': // msft
              if (tmpstr == NULL || sscanf(tmpstr, "%d", &device_number) != 1)
                device_number = 0;
              else
                if (tmpstr != NULL && tmpstr == argv[arg + 1])
                  arg++;
              break;
       /* case 'd': is with case 'o': below because they share code */
          case 'h': /* human-readable checkpoint format */
            if (strncmp(program_name, "fftlucas", 8))
            {
              fprintf(stderr, "%s: only fftlucas supports human readable checkpoints\n%s", program_name, Options);
//              exit(1);
            }
            chkpnt_fmt = CHKPNT_ASCII;
            break;
          case 'o': case 'd': /* 'o'utput and 'd'uplicate output files */
            i = argv[arg][1];
            if (tmpstr == NULL)
            {
              fprintf(stderr, "%s: no file name after -%c\n%s", program_name, i, Options);
              exit(1);
            }
            if (!strcmp(tmpstr, "stdout") || /* rest of condition is "-o-", "-o -", "-o -anythingelse" */
                (i == 'o' && (!strcmp(tmpstr, "-") || (argv[arg][2] == '\0' && argv[arg + 1][0] == '-'))))
              tmpfp = stdout; /* that 'rest' of the condition avoids "-ofile -n", "-o-file", etc */
            else
              if (!strcmp(tmpstr, "stderr") ||
                  (i == 'd' && (!strcmp(tmpstr, "-") || (argv[arg][2] == '\0' && argv[arg + 1][0] == '-'))))
                tmpfp = stderr;
              else
                if ((tmpfp = fopen(tmpstr, "a")) == NULL)
                {
                  fprintf(stderr, "%s: could not open '%s' for append\n", program_name, tmpstr);
                  perror(tmpstr);
                  exit(errno);
                }
            if (tmpstr == argv[arg + 1] && (tmpstr[0] != '-' || tmpstr[1] == '\0'))
              arg++;
            if (i == 'o')
            { /* normal output */
              if (*outfp != NULL && *outfp != stdout && *outfp != stderr)
                fclose(*outfp);
              *outfp = tmpfp;
              if (*dupfp == *outfp)
                *dupfp = NULL; /* no point in duplicating the output to the same descriptor */
            }
            else
            { /* duplicate output, usually to get local copy when the normal output is via TCP/IP */
              if (*dupfp != NULL && *dupfp != stdout && *dupfp != stderr)
                fclose(*dupfp);
              if (tmpfp == *outfp)
                tmpfp = NULL;
              *dupfp = tmpfp;
            }
            tmpfp = NULL;
            break;
          case 'r':
            if (tmpstr == NULL || !isascii(tmpstr[0]) || !isdigit(tmpstr[0]) ||
                (tmpstr == argv[arg + 1] &&
                 (argv[arg + 2] == NULL || !isascii(argv[arg + 2][0]) || !isdigit(argv[arg + 2][0]))) ||
                (tmpstr != argv[arg + 1] &&
                 (argv[arg + 1] == NULL || !isascii(argv[arg + 1][0]) || !isdigit(argv[arg + 1][0]))))
            {
              fprintf(stderr, "%s: -r requires two numeric arguments\n%s", program_name, Options);
              exit(1);
            }
            sscanf(tmpstr, SCANF_FMT_UL, &start_range);
            if (tmpstr == argv[++arg])
              sscanf(argv[++arg], SCANF_FMT_UL, &stop_range);
            else
              sscanf(argv[arg], SCANF_FMT_UL, &stop_range);
            if (start_range < (UL)31 || stop_range > (~(UL)0)/2 || start_range > stop_range)
             {
              fprintf(stderr, "%s: -r start too small, stop too large, or start > stop\n", program_name);
              exit(1);
             }
            break;
          case 't':
            print_times = 1;
            break;
          case 'v':
            fprintf(stderr, "%s version information:\n %s\n", program_name, version_info);
            exit(0);
          default:
            fprintf(stderr, "%s: unknown option: '%s'\n%s", program_name, argv[arg], Options);
            exit(1);
         }
        arg++;
        break;
      default: /* anything except a leading '-' */
        tmpstr = argv[arg];
        char buffer[CHKPNT_FILENAME_LENGTH];
        tmpfp = fopen(tmpstr, "rb");
        if (tmpfp == NULL)
        {
          sprintf(buffer, "c%s", tmpstr);
          tmpfp = fopen(buffer,"rb");
          if (tmpfp == NULL)
          {
            sprintf(buffer, "t%s", tmpstr);
            tmpfp = fopen(buffer, "rb");
            if (tmpfp != NULL)
              fprintf(stdout, "%s: Resuming from backup checkpoint file %s\n", program_name, buffer);
            else
              fprintf(stderr, "%s: Could not find a checkpoint file to resume from\n", program_name);
          }
          else
            fprintf(stdout, "%s: Resuming from checkpoint file %s\n", program_name, buffer);
         }
         else
           fprintf(stdout, "%s: Resuming from checkpoint file %s\n", program_name, tmpstr);
 
         if (tmpfp != NULL)
         {
          *infp = tmpfp;
          first_time_here = 0;
          tmpfp = NULL;
         }
        else
          if (sscanf(argv[arg], SCANF_FMT_UL, &tmp) == 1)
           {
            *q = tmp;
            *infp = NULL;
            first_time_here = 0;
           }
          else
           {
            fprintf(stderr, "%s: tried to read '%s' as a check point file or list of exponents\n%s",
                          program_name, argv[arg], Options);
            perror(program_name);
            exit(errno == 0 ? 1 : errno);
           }
        break;
     }
  }
  if (*infp == NULL && first_time_here) *infp = stdin;
  if (*outfp == NULL) *outfp = stdout;
  first_time_here = false;

  if (*infp != NULL && fscanf(*infp, " trial%c", str) == 1 && str[0] == 's')
  {
    /* gobble the digits, including any leading whitespace and escaped linefeeds and returns */
    while (((i = fgetc(*infp)) != EOF && isascii(i) && isspace(i)) || (i == '\\' && ((i = fgetc(*infp)) == '\n' || i == '\r')));
    while ((i != EOF && isascii(i) && isdigit(i)) || (i == '\\' && ((i = fgetc(*infp)) == '\n' || i == '\r')))
      i = fgetc(*infp);
  }   

  if (*infp == NULL || fscanf(*infp, SCANF_FMT_UL, q) == 1)
  { /* if there's no input file, *q is from the command line */
    if (*q < (UL)31)
    {
      fprintf(stderr, "%s: exponent " PRINTF_FMT_UL " is too small to test\n", program_name, *q);
      return(1);
    }
    if (*infp != NULL)
      switch (i = fgetc(*infp))
      {
        default: /* unknown format; presume it's just a list of exponents */
          ungetc(i, *infp);
          /* fall thru */
        case '\n': /* list of exponents */
        case ' ':
          *M = 'U';
          break;
      }
    if (n != NULL)
    { /* must be LL testing, not factoring */
      if (*q < start_range || *q > stop_range)
        goto get_next_q;
      *M = 'U';
      if (last >= *q) /* reset n in case it's too big for this new exponent */
      {
          *n = 8L;
      }
      sprintf(chkpnt_cfn, "c" PRINTF_FMT_UL, *q);
      sprintf(chkpnt_tfn, "t" PRINTF_FMT_UL, *q);
      /*!!see if the checkpoint files already exist and continue from them? */
      /*!!can't easily use *infp to do it since that would lose track of where we are in the current input */
      /*!!would help to have new read_chkpoint() and save_chkpoint() functions */
      /*!!new good_chkpoint() function that returns flag: is checkpoint data just read consistent? */
      return(3);
    }
    return(1); /* incorrect call since we must be doing LL tests */
  }

  if (fscanf(*infp, "M(" SCANF_FMT_UL " )%c", q, M) == 2)
  {
    if (*q < (UL)31 || presence("UHIEeoqGCDdPJcLMlm", *M) == NULL)
    {
      fprintf(stderr, "%s: invalid data (exponent or letter 'code')\n", program_name);
      return(1); /* different format or error */
    }
    if (n != NULL)
    {
      while ((i = fgetc(*infp)) != EOF && i != '\n')
        ; /* rest of line is useless to LL testers */ /*!!unless to make sure it's factored far enough */
      if (*q < start_range || *q > stop_range)
        goto get_next_q;
      if (last >= *q) /* reset n in case it's too big for this new exponent */
      {
          *n = 8L;
      }
      sprintf(chkpnt_cfn, "c" PRINTF_FMT_UL, *q);
      sprintf(chkpnt_tfn, "t" PRINTF_FMT_UL, *q);
      return(3);
    }
    return(1); /* incorrect call */
  }

  if (feof(*infp)) /* this was already checked except in the case of whitespace */
    goto get_next_q; /*  - including a solitary return or LF - in front of EOF */
  /* everything else is for Lucas-Lehmer testers */
  if (x == NULL) /* no place to store pointer to FFT array from checkpoint */
    return(1); /* then the caller goofed; let them know */
  /* check for printable/human readable checkpoint header; partial will become 7 if so */
  partial = fscanf(*infp, SCANF_FMT_ONE, q, n, j, err, &sBL, &sD, &RHM);
  if (partial < 0)
  {
    if (ferror(*infp)) return 1;
    else goto get_next_q; /* would be return(0) except that there might be more command line to parse */
  }

  if (partial == 0 && !ferror(*infp))
  { /* check for binary checkpoint format */
    int r = read_check_point(*infp, q, n, j, err, x);
    if(r)
    {
      *infp = NULL;
      return r;
    }
  }

  if (partial != 7 || *n < ((double)(*q))/ROUNDED_HALF_MANTISSA + 4 || *q < (UL)31 || *q % 2 == 0 ||
      *j >= *q || *j < 1 || *err > 0.49 || (*err > 0.3 && *q < (UL)10000) || (sizeof(BIG_LONG) != sBL &&
      sBL != 4 && sBL != 8) || sizeof(double) != sD || ROUNDED_HALF_MANTISSA != RHM)
  {
    fprintf(stderr, "%s: inconsistent RI header (initial line)\n", program_name);
    return(1);
  }

  for (tmp2 = *n; (tmp2 & 1) == 0; tmp2 >>= 1);

  if (tmp2 != 1 && tmp2 != 3 && tmp2 != 5 && tmp2 != 7)
  {
    fprintf(stderr, "%s: inconsistent RI input; FFT run length not a power of 2\n", program_name);
    return(1);
  }
  while (*n < ((double)(*q))/ROUNDED_HALF_MANTISSA + 4.0)
    *n <<= 1;
  if (*x != NULL)
    free((char *)*x);

  if (((*x) = (double *)calloc((size_t)(*n + *n), sizeof(double))) == NULL)
  {
    perror(program_name);
    fprintf(stderr, "%s: cannot get memory for FFT array\n", program_name);
    return(1);
  }
  for (tmp2 = 0; tmp2 < *n; tmp2++)
    if (fscanf(*infp, SCANF_FMT_TWO, x[0] + tmp2) != 1)
    {
      perror(program_name);
      fprintf(stderr, "%s: error reading RI input (middle lines)\n", program_name);
      return(1);
  }

  if (fscanf(*infp, " RI " PRINTF_FMT_UL " done", &tmp2) != 1)
  {
    perror(program_name);
    fprintf(stderr, "%s: error reading RI input (last line)\n", program_name);
    return(1);
  }
  if (*q != tmp2)
  {
    fprintf(stderr, "%s: inconsistent RI input (last line)\n", program_name);
    return(1);
  }

  sprintf(chkpnt_cfn, "c" PRINTF_FMT_UL, *q);
  sprintf(chkpnt_tfn, "t" PRINTF_FMT_UL, *q);
  fclose (*infp);

  return 2; /*!!could get_next_q here also, but skip an incomplete checkpoint?? */
}


/* print what we've done so far to a (local) checkpoint file */
int check_point(UL q, UL n, UL j, double err, double *x)
{
  FILE *chkpnt_fp;
  UL k = MAGIC_NUMBER; /* magic number ala the file(1) command; note that it's byte order dependent, */
                        /* which we want so that we don't try to use a save file from a different */
                        /* architecture or byte order */
  size_t retval; /* return value from fwrite(); Linux man page and libc.a do not agree on correct value */

  /* refuse to write checkpoint files for small exponents */
  if (q < 1000) return 1 ;

  if ((chkpnt_fp = fopen(chkpnt_cfn, "rb")) != NULL)
  {
    fclose(chkpnt_fp);
    UNLINK(chkpnt_tfn);
    rename(chkpnt_cfn, chkpnt_tfn);
  }
  if ((chkpnt_fp = fopen(chkpnt_cfn, "wb")) == NULL)
  {
    perror(program_name);
    fprintf(stderr, "%s: could not open '%s' for writing\n", program_name, chkpnt_cfn);
    return(0);
  }
  switch (chkpnt_fmt)
  {
    case CHKPNT_BINARY:
      /* 'tag' the file as coming from this program, computer architecture, etc. */
      /* Note that k was initialized to MAGIC_NUMBER as fwrite() can't write a #define. */
      if ((retval = fwrite(&k, sizeof(k), 1UL, chkpnt_fp)) != sizeof(k)
       && retval != 1UL)
      {
        perror(program_name);
        fprintf(stderr, "%s: cannot write checkpoint info (k = " PRINTF_FMT_UL ")\n", program_name, k);
        fprintf(stderr, "%s: fwrite returned %d\n", program_name, (int)retval);
        return(0);
      }
      if ((retval = fwrite(&q, sizeof(q), 1, chkpnt_fp)) != sizeof(q) && retval != 1UL)
      {
        perror(program_name);
        fprintf(stderr, "%s: cannot write checkpoint info (q = " PRINTF_FMT_UL ")\n", program_name, q);
        fprintf(stderr, "%s: fwrite returned %d\n", program_name, (int)retval);
        return(0);
      }
      if ((retval = fwrite(&n, sizeof(n), 1, chkpnt_fp)) != sizeof(n) && retval != 1UL)
      {
        perror(program_name);
        fprintf(stderr, "%s: cannot write checkpoint info (n = " PRINTF_FMT_UL ", q = " PRINTF_FMT_UL ")\n",
                program_name, n, q);
        fprintf(stderr, "%s: fwrite returned %d\n", program_name, (int)retval);
        return(0);
      }
      if ((retval = fwrite(&j, sizeof(j), 1, chkpnt_fp)) != sizeof(j) && retval != 1UL)
      {
        perror(program_name);
        fprintf(stderr, "%s: cannot write checkpoint info (j = " PRINTF_FMT_UL ", q = " PRINTF_FMT_UL ")\n",
                program_name, j, q);
        fprintf(stderr, "%s: fwrite returned %d\n", program_name, (int)retval);
        return(0);
      }
      if ((retval = fwrite(&err, sizeof(err), 1, chkpnt_fp)) != sizeof(err) && retval != 1UL)
      {
        perror(program_name);
        fprintf(stderr, "%s: cannot write checkpoint info (err, q = " PRINTF_FMT_UL ")\n", program_name, q);
        fprintf(stderr, "%s: fwrite returned %d\n", program_name, (int)retval);
        return(0);
      }
      for (k = 0; k < n; k++)
        if ((retval = fwrite(&(x[k]), sizeof(x[0]), 1, chkpnt_fp)) != sizeof(x[0]) && retval != 1UL)
        {
          perror(program_name);
          fprintf(stderr, "%s: cannot write checkpoint info (x, q = " PRINTF_FMT_UL ")\n", program_name, q);
          fprintf(stderr, "%s: fwrite returned %d\n", program_name, (int)retval);
          return(0);
        }
      break;
    case CHKPNT_ASCII:
    default: /* some compilers think sizeof() returns an int and others think it returns a long */
      if (fprintf(chkpnt_fp, PRINTF_FMT_ONE, q, n, j, err, (long)sizeof(BIG_LONG), (long)sizeof(double),
                  ROUNDED_HALF_MANTISSA) <= 0)
      {
        perror(program_name);
        fprintf(stderr, "%s: cannot write checkpoint info (RI header, q = " PRINTF_FMT_UL ")\n",
                program_name, q);
        return(0);
      }
      for (k = 0; k < n; k++)
        if (fprintf(chkpnt_fp, PRINTF_FMT_TWO, x[k]) <= 0)
        {
          perror(program_name);
          fprintf(stderr, "%s: cannot write checkpoint info (RI mid, q = " PRINTF_FMT_UL ")\n",
                  program_name, q);
          return(0);
        }
      if (fprintf(chkpnt_fp, "RI " PRINTF_FMT_UL " done\n", q) <= 0)
      {
        perror(program_name);
        fprintf(stderr, "%s: cannot write checkpoint info (RI tail, q = " PRINTF_FMT_UL ")\n",
                program_name, q);
        return(0);
      }
      break;
  }
  fflush(chkpnt_fp);
#if !defined(macintosh) && !defined(_MSC_VER)
  fsync(fileno(chkpnt_fp));
#endif
  return(!fclose(chkpnt_fp));
}

static void note_archive_failure(FILE *fp, int error)
{
  static char cnt = 0;

  if (++cnt > 3)
    return; /* already warned in all possible output files, but warn again after it overflows */
  if (fp != NULL && !ferror(fp))
  {
    fprintf(fp, "%s: previous line may not have been appended to %s correctly\n", program_name, archive_name());
    fprintf(fp, "\terrno %d: %s\n", error, strerror(error));
  }
}

static void close_archive(FILE *archfp, FILE *outfp, FILE *dupfp)
{
  int last_errno;

  if (fclose(archfp) != 0)
  {
    note_archive_failure(outfp, last_errno = errno);
    note_archive_failure(dupfp, last_errno);
    if (outfp != stderr && dupfp != stderr)
      note_archive_failure(stderr, last_errno);
    return;
  }
  UNLINK(chkpnt_cfn);
  UNLINK(chkpnt_tfn);
}

void printbits(double *x, UL q, UL n, UL totalbits, UL b, UL c, double hi, double lo,
     const char *version_info, FILE *outfp, FILE *dupfp, int iterations, int current_iteration, bool archive)
{
  FILE *archfp = NULL;

  if(archive) archfp = open_archive();

  if (version_info == NULL || version_info[0] == '\0')
  {
    fprintf(stderr, "%s: printbits() called with no version info; bug\n", program_name);
    if (outfp != NULL && outfp != stderr)
    {
      fprintf(outfp, "%s: printbits() called with no version info; bug\n", program_name);
    }
    if (dupfp != NULL && dupfp != stderr)
    {
      fprintf(dupfp, "%s: printbits() called with no version info; bug\n", program_name);
    }
  }

  if(!archive) fprintf(outfp,"Iteration %d ", current_iteration);
  if (is_zero(x, n, 0, 0))
  {
    if (outfp != NULL && !ferror(outfp)) fprintf(outfp,  "M( " PRINTF_FMT_UL " )P, n = " PRINTF_FMT_UL ", %s\n", q, n, version_info);
    if (dupfp != NULL && !ferror(dupfp)) fprintf(dupfp,  "M( " PRINTF_FMT_UL " )P, n = " PRINTF_FMT_UL ", %s\n", q, n, version_info);
    if (archfp != NULL)                  fprintf(archfp, "M( " PRINTF_FMT_UL " )P, n = " PRINTF_FMT_UL ", %s\n", q, n, version_info);

    if (archfp != NULL) close_archive(archfp, outfp, dupfp);

    return;
  }

  if (outfp != NULL && !ferror(outfp)) fprintf(outfp,  "M( " PRINTF_FMT_UL " )C", q);
  if (dupfp != NULL && !ferror(dupfp)) fprintf(dupfp,  "M( " PRINTF_FMT_UL " )C", q);
  if (archfp != NULL)                  fprintf(archfp, "M( " PRINTF_FMT_UL " )C", q);

  if (totalbits < 1)
  {
    if (outfp != NULL && !ferror(outfp)) fprintf(outfp,  ", %s\n", version_info);
    if (dupfp != NULL && !ferror(dupfp)) fprintf(dupfp,  ", %s\n", version_info);
    if (archfp != NULL)                  fprintf(archfp, ", %s\n", version_info);

    if (archfp != NULL) close_archive(archfp, outfp, dupfp);

    return;
  }
  balancedtostdrep(x, n, b, c, hi, lo, 0, 0);
  static UL *hex = NULL;
  static UL prior_hex = 0;
  if (hex != NULL && totalbits/BITS_IN_LONG + 1 > prior_hex)
  {
    free(hex);
    hex = NULL;
    prior_hex = totalbits/BITS_IN_LONG + 1;
  }
  if (hex == NULL && (hex = (UL *)calloc(totalbits/BITS_IN_LONG + 1, sizeof(UL))) == NULL)
  {
    perror(program_name);
    fprintf(stderr, "%s: cannot get memory for residue bits; calloc() errno %d\n", program_name, errno);
    exit(errno);
  }

  int i = 0;
  int j = 0;
  do
  {
    UL k = (long)(ceil((double)q*(j + 1)/n) - ceil((double)q*j/n));
    if (k > totalbits) k = totalbits;
    totalbits -= k;
    int word = (long)x[j];
    for (j++; k > 0; k--, i++)
    {
      if (i % BITS_IN_LONG == 0) hex[i/BITS_IN_LONG] = 0L;
      hex[i/BITS_IN_LONG] |= ((word & 0x1) << (i % BITS_IN_LONG));
      word >>= 1;
    }
  } while(totalbits > 0);

  if (outfp != NULL && !ferror(outfp)) fputs(", 0x", outfp);
  if (dupfp != NULL && !ferror(dupfp)) fputs(", 0x", dupfp);
  if (archfp != NULL)                  fputs(", 0x", archfp);

  static char bits_fmt[16] = "\0"; /* "%%0%ulx" -> "%08lx" or "%016lx" depending on sizeof(UL) */

  if (bits_fmt[0] != '%')
    sprintf(bits_fmt, "%%0%lu%s", (UL)(BITS_IN_LONG/4), "lx"); /* 4 bits per hex 'digit' */

  for (j = (i - 1)/BITS_IN_LONG; j >= 0; j--)
  {
    if (outfp != NULL && !ferror(outfp)) fprintf(outfp, bits_fmt, hex[j]);
    if (dupfp != NULL && !ferror(dupfp)) fprintf(dupfp, bits_fmt, hex[j]);
    if (archfp != NULL)                  fprintf(archfp, bits_fmt, hex[j]);
  }

  if (outfp != NULL && !ferror(outfp))
  {
    fprintf(outfp, ", n = " PRINTF_FMT_UL ", %s", n, version_info);
    if(!archive && print_times)
    {
      printf(" (");
      print_time(iterations, current_iteration, q);
      printf(")");
    }
    fprintf(outfp, "\n");
    fflush(outfp);
  }
  if (dupfp != NULL && !ferror(dupfp))
    fprintf(dupfp, ", n = " PRINTF_FMT_UL ", %s\n", n, version_info);
  if (archfp != NULL)
  {
    fprintf(archfp, ", n = " PRINTF_FMT_UL ", %s\n", n, version_info);
    close_archive(archfp, outfp, dupfp);
  }
  return;
}
