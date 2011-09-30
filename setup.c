static const char RCSsetup_c[] = "$Id: setup.c,v 8.1 2007/06/23 22:33:35 wedgingt Exp $";

#include "setup.h"

static const char RCSsetup_h[] = DRCSsetup_h;

#ifndef vms
volatile char shouldTerminate = 0; /* Flag: have we gotten a SIGTERM ? */
#endif

/* only used to make the process's priority as low as possible */
//#if !defined(macintosh) && !defined(_MSC_VER)
//extern int nice (int priority_increment);
//#endif

#ifdef __sgi__
#include <limits.h>
#include <sys/types.h>
#include <sys/prctl.h>
#include <sys/schedctl.h>
#endif

/* same as the BSD index() and the SYSV strchr() */
const char *presence(const char *string, int ch)
{
  if (string == NULL)
    return(NULL);
  while (*string != '\0')
    if (*string == ch)
      return(string);
    else
      ++string;
  return(NULL);
}

void setup (void)
 {
  setlinebuf(stdout);
//  setlinebuf(stderr);
#ifdef __sgi__
  (void)schedctl(NDPRI, 0, NDPLOMIN);
#endif
#if !defined(macintosh) && !defined(_MSC_VER)
//  (void)nice(40);
#endif
#ifdef SIGTERM
  (void)signal(SIGTERM, term_handler);
#endif
#ifdef SIGINT
  (void)signal(SIGINT, term_handler);
#endif
#ifdef SIGHUP
  (void)signal(SIGHUP, term_handler);
#endif
#ifdef SIGPIPE
  (void)signal(SIGPIPE, SIG_IGN);
#endif
#ifdef ANCIENT_linux
  /*!!Should not be needed any more, as implied by the new macro protecting it.  Old comments: */
  /*!!really ought to examine the signal and let it happen if it's an overflow, NaN, etc.,*/
  /*!! but it's almost always an underflow or round off warning rather than a problem like overflow or NaN */
  (void)signal(SIGFPE, SIG_IGN);
#endif
  return;
 }

#if !defined(vms)
/*ARGSUSED*/ /* they aren't but that's because only some OS's want an argument here */
# if defined(linux) || defined(__ultrix) || defined(_AIX) || defined(__hpux) || defined(macintosh) || defined(__APPLE__) || defined(_MSC_VER)
handler term_handler(int)
# else
handler term_handler()
# endif
 {
# ifdef SIGTERM
  (void)signal(SIGTERM, SIG_IGN);
# endif
# ifdef SIGINT
  (void)signal(SIGINT, SIG_IGN);
# endif
# ifdef SIGHUP
  (void)signal(SIGHUP, SIG_IGN);
# endif
  shouldTerminate = 1;
  return_handler;
}

# ifndef pc7300
void clientexit(const char *msg)
{
  (void)fprintf(stderr, "%s: ", program_name);
  (void)fflush(stderr);
  perror(msg);
  exit(errno);
}

#  ifdef NO_HERROR
#   ifndef NO_ADDRESS
#    define NO_ADDRESS 3
#   endif
void herror(const char *msg)
{
  static const char *(herrs[]) = {
    "unknown error",
    "no such host",
    "server failure or host not known",
    "no address for host"
  };

  (void)fprintf(stderr, "%s: %s (code %d)\n", msg, h_errno <= NO_ADDRESS ? herrs[h_errno] : herrs[0],
        h_errno);
}
#  endif
# endif
#endif
