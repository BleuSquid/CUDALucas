#ifndef DRCSsetup_h
#define DRCSsetup_h "$Id: setup.h,v 8.1 2007/06/23 22:33:35 wedgingt Exp $"

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

extern const char *program_name;
extern const char program_revision[];

#define BIG_LONG long int
#define BIG_LONG_FMT "%ld"
#define SIZEOF_UL 4
  /* the previous would not be needed if the preprocessor understood sizeof() */
#define ROUNDED_HALF_MANTISSA 32L
/* truncate an IEEE double by adding and then subtracting 3.0*2^51 to it */
#define TRUNC_A (3.0*0x4000000*0x2000000)
#define TRUNC_B (12.0*0x2000000*0x1000000)
#define PRINTF_FMT_ONE "RI q %lu n %lu j %lu err %f sBL %ld sD %ld RHM %ld\n"
#define PRINTF_FMT_TWO "RI %.1f\n"
#define PRINTF_FMT_UL "%lu"
#define SCANF_FMT_UL " %lu"
#define SCANF_FMT_ONE " RI q %lu n %lu j %lu err %lf sBL %d sD %d RHM %ld\n"
#define SCANF_FMT_TWO " RI %lf\n"

# include <stdlib.h>

/* This has to come before the #define's of exp() et al on some Linux machines */
# include <math.h>

# if defined(__GNUC__) && !defined(NeXT)
#  if (__GNUC__ == 2 && __GNUC_MINOR__ >= 6) || __GNUC__ > 2
#   if defined(__i386__) && !defined(__CYGWIN__)
#    ifndef DO_NOT_USE_LONG_LONG_EXPONENT
#     undef BIG_LONG
#     define BIG_LONG long long int
#     undef SIZEOF_UL
#     define SIZEOF_UL 8
#    endif
#   endif
#  endif
# endif

#if defined(__GNUC__)
# define ATTRIBUTE_PURE		__attribute__((pure))
# define ATTRIBUTE_CONST	__attribute__((const))
# define ATTRIBUTE_NORETURN	__attribute__((noreturn))
#else
# define ATTRIBUTE_PURE
# define ATTRIBUTE_CONST
# define ATTRIBUTE_NORETURN
#endif

# include <stdio.h>
# include <ctype.h>

# include <assert.h>

# if !defined(__amiga__) && !defined(macintosh) && !defined(__FreeBSD__) && !defined(__APPLE__)
#  include <malloc.h>
# endif

# ifdef rint
#  undef rint
# endif
# define rint(x) ((double)((BIG_LONG)(x + 0.5)))

#  include <string.h>

# ifndef pc7300

#  ifdef _MSC_VER
#   define bzero(str, len) (void)memset(str, '\0', len)
#  endif

# else
#  include <memory.h>

#  define bzero(str, len) (void)memset(str, '\0', len)

/* CAREFUL: this macro only works when newfd < smallest fd not open */
/*  We only need it to replace 0 (stdin) and 1 (stdout) presently */
/*  It could be done properly with pair of while loops and an array to */
/*   remember the fd's we opened that we didn't want */
#  define dup2(fd, newfd) ((void)close(newfd), (void)dup(fd))

#  define setlinebuf(fp) setbuf(fp, NULL)

# endif

# include <signal.h>
# include <errno.h>

#ifdef _IOLBF
# define setlinebuf(stream) setvbuf(stream, NULL, _IOLBF, BUFSIZ)
// Why do we need this here? We've already included stdio.h.
//extern int setvbuf (FILE *stream, char *buffer, int type, size_t size);
#endif

typedef unsigned short US;

typedef unsigned long int UL;

typedef void handler;
#define return_handler return

extern volatile char shouldTerminate;   /* Flag: have we gotten a SIGTERM ? */

/* NBBY is Number of Bits per BYte, usually 8 */
#ifndef NBBY
# define NBBY 8
#endif

#define BITS_IN_LONG (NBBY*sizeof(UL))

#define DIGITS_IN_LONG ((BITS_IN_LONG + 1)*10/31)
/* 32 bits fit in 9.6+ digits */

#define CHKPNT_FILENAME_LENGTH (DIGITS_IN_LONG + 2)
/* + 2 for c (or t or ...) and \0 */

ATTRIBUTE_PURE extern const char *presence (const char *string, int ch); /* same as BSD index() and SysV strchr() */
extern void setup (void);

# if defined(linux) || defined(__ultrix) || defined(_AIX) || defined(__hpux) || defined(macintosh) || defined(_MSC_VER) || defined(__APPLE__) || defined(_MSC_VER)
extern handler term_handler (int huh);
# else
extern handler term_handler (void);
# endif

#endif /* ifndef DRCSsetup_h */
