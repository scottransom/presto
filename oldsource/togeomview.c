
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/signal.h>
#include <sys/file.h>
#ifdef AIX
#define _BSD  1			/* Get FNDELAY from fcntl.h */
#endif
#include <fcntl.h>
#include <sys/errno.h>
#include <sys/socket.h>
#include <sys/un.h>

#ifndef NeXT
#include <unistd.h>
#endif

#ifndef O_NONBLOCK
#define O_NONBLOCK O_NDELAY
#endif
#ifndef FNDELAY
#define FNDELAY O_NDELAY
#endif
#ifndef FNONBLK
#define FNONBLK FNDELAY
#endif

extern int errno;

char *todir = "/tmp/geomview";
char *toname = "OOGL";
char *prog;

char timedout[] = "Couldn't start geomview\n";

static void interrupt(int sig)
{
  write(2, timedout, sizeof(timedout));
  exit(1);
}

void usage(void);

/*
 * We've determined that no copy of geomview is running, so start one.
 */
void start_gv(char **progtorun, char *toname, int asgeom, int usesock)
{
  char *args[1024];
  int i;
  static char Mprefix[] = "/tmp/geomview/";
  static char Mhow[] = "-Mcp";

  signal(SIGALRM, interrupt);
  for (i = 0; progtorun[i] != NULL; i++)
    args[i] = progtorun[i];

  Mhow[2] = asgeom ? 'g' : 'c';
  Mhow[3] = usesock ? 's' : 'p';
  args[i++] = Mhow;
  args[i++] = strncmp(toname, Mprefix, strlen(Mprefix))
      ? toname : toname + strlen(Mprefix);

  if (fork() == 0) {
    close(0);
#if defined(NeXT) || defined(BSD)
    setpgrp(0, getpid());
#else				/* SysV style */
    setpgrp();
#endif
    execvp(progtorun[0], &args[0]);
    execvp("geomview", &args[0]);
    perror("Couldn't exec geomview");
    kill(getppid(), SIGALRM);
    _exit(1);
  }
}

/*
 * Pipe fitting for linking Mathematica to geomview.
 * Starts geomview if not already running.
 */
int main(int argc, char *argv[])
{
  int usesock;
  int n, fd = -1;
  int asgeom = 0;
  char pipename[BUFSIZ];
  static char *geomview[] =
  {"geomview", NULL};
  char **progtorun = geomview;
  char *tail;
  struct sockaddr_un un;
  struct stat st;

  prog = argv[0];
  tail = strrchr(prog, '/');
  if (tail)
    tail++;
  else
    tail = argv[0];
  if (tail[0] == 'g') {
    asgeom = 1;
  }
#ifdef NeXT
  usesock = 1;
#else
  usesock = 0;
#endif

  while (argc > 1 && argv[1][0] == '-') {
    for (tail = argv[1]; *++tail;) {
      switch (*tail) {
      case 'M':
	asgeom = 1;
	break;
      case 'g':
	asgeom = 1;
	break;
      case 'c':
	asgeom = 0;
	break;
      case 'p':
	usesock = 0;
	break;
      case 's':
	usesock = 1;
	break;
      default:
	usage();
      }
    }
    argc--, argv++;
  }

  if (argc > 1)
    toname = argv[1];
  if (argc > 2)
    progtorun = &argv[2];

  if (toname[0] == '/') {
    tail = strrchr(toname, '/');
    *tail = '\0';
    todir = toname;
    toname = tail + 1;
  }
  if (access(todir, W_OK) < 0) {
    mkdir(todir, 0777);
    chmod(todir, 0777);
  }
  sprintf(pipename, "%s/%s", todir, toname);

  if (stat(pipename, &st) >= 0 &&
      (st.st_mode & S_IFMT) != (usesock ? S_IFSOCK : S_IFIFO))
    unlink(pipename);

  if (usesock) {
    strncpy(un.sun_path, pipename, sizeof(un.sun_path) - 1);
    un.sun_family = AF_UNIX;
    fd = socket(AF_UNIX, SOCK_STREAM, 0);
    if (connect(fd, (struct sockaddr *) &un, sizeof(un)) < 0) {
      if (errno != ECONNREFUSED && errno != ENOENT) {
	fprintf(stderr, "togeomview: Can't connect to ");
	perror(pipename);
	exit(1);
      }
      start_gv(progtorun, pipename, asgeom, usesock);
      for (n = 0; connect(fd, (struct sockaddr *) &un, sizeof(un)) < 0; n++) {
	if (n == 15)
	  interrupt(0);
	sleep(1);
      }
    }
  } else {
    /* Use named pipe */
#ifdef S_IFIFO

    if (access(pipename, 0) < 0) {
      mknod(pipename, S_IFIFO, 0);
      chmod(pipename, 0666);
    }
    fd = open(pipename, O_WRONLY | O_NONBLOCK);
    if (fd >= 0) {
      fcntl(fd, F_SETFL, fcntl(fd, F_GETFL, 0) & ~(FNDELAY | FNONBLK));
    } else if (errno == ENXIO) {
      start_gv(progtorun, pipename, asgeom, usesock);
      alarm(60);
      fd = open(pipename, O_WRONLY);
      alarm(0);
    }
#endif
  }
  if (fd < 0) {
    fprintf(stderr, "Can't open pipe to geomview: ");
    perror(pipename);
    exit(1);
  }
  while ((n = read(0, pipename, sizeof(pipename))) > 0) {
    if (write(fd, pipename, n) < n) {
      perror("Error writing to geomview");
      exit(1);
    }
  }
  exit(0);
}

void usage(void)
{
  char buf[BUFSIZ];
  setbuf(stderr, buf);

  fprintf(stderr, "Usage: %s [-g] [pipename  [ geomview ... args ] ]\n\
Sends lisp-style commands or (with \"-g\") OOGL geometry data to geomview,\n\
starting a copy if none is yet running.  Uses \"pipename\" as the connection\n\
name; a file by that name is created in the directory \"/tmp/geomview\".\n\
Default pipename is \"OOGL\".  If \"geomview ... args\" are present,\n\
invokes that (if need be) rather than geomview itself.\n\
Examples:\n\
	echo '(geometry fred < dodec.off)' | togeomview  sam\n\
	togeomview -g <dodec.off\n\
	cat my_geomview_script | togeomview bob  geomview -wpos 300x300@500,500\n",
	  prog);
  exit(1);
}
