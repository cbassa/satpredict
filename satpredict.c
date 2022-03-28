#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include "sgdp4h.h"

#define LIM 256
#define NMAX 16
#define D2R M_PI/180.0
#define R2D 180.0/M_PI
#define XKMPER 6378.135 // Earth radius in km
#define XKMPAU 149597879.691 // AU in km
#define FLAT (1.0/298.257)

long Isat=0;
long Isatsel=0;
extern double SGDP4_jd0;

struct map {
  double lat,lng;
  float alt;
  char observer[32];
  int site_id;
} m;
struct point {
  double mjd;
  xyz_t obspos,sunpos;
  double zeta,z,theta;
} *p;
double modulo(double,double);
void obspos_xyz(double,xyz_t *,xyz_t *);
void sunpos_xyz(double,xyz_t *);
double gmst(double);
double dgmst(double);
void precession_angles(double mjd0,double mjd,double *zeta,double *z,double *theta);
double nfd2mjd(char *date);
double date2mjd(int year,int month,double day);

// Compute Julian Day from Date
double date2mjd(int year,int month,double day)
{
  int a,b;
  double jd;

  if (month<3) {
    year--;
    month+=12;
  }

  a=floor(year/100.);
  b=2.-a+floor(a/4.);

  if (year<1582) b=0;
  if (year==1582 && month<10) b=0;
  if (year==1582 && month==10 && day<=4) b=0;

  jd=floor(365.25*(year+4716))+floor(30.6001*(month+1))+day+b-1524.5;

  return jd-2400000.5;
}

// nfd2mjd
double nfd2mjd(char *date)
{
  int year,month,day,hour,min;
  float sec;
  double mjd,dday;

  sscanf(date,"%04d-%02d-%02dT%02d:%02d:%f",&year,&month,&day,&hour,&min,&sec);
  dday=day+hour/24.0+min/1440.0+sec/86400.0;
  mjd=date2mjd(year,month,dday);

  return mjd;
}

// Greenwich Mean Sidereal Time
double gmst(double mjd)
{
  double t,gmst;

  t=(mjd-51544.5)/36525.0;

  gmst=modulo(280.46061837+360.98564736629*(mjd-51544.5)+t*t*(0.000387933-t/38710000),360.0);

  return gmst;
}

// Greenwich Mean Sidereal Time
double dgmst(double mjd)
{
  double t,dgmst;

  t=(mjd-51544.5)/36525.0;

  dgmst=360.98564736629+t*(0.000387933-t/38710000);

  return dgmst;
}

// Return x modulo y [0,y)
double modulo(double x,double y)
{
  x=fmod(x,y);
  if (x<0.0) x+=y;

  return x;
}

// Observer position
void obspos_xyz(double mjd,xyz_t *pos,xyz_t *vel)
{
  double ff,gc,gs,theta,s,dtheta;

  s=sin(m.lat*D2R);
  ff=sqrt(1.0-FLAT*(2.0-FLAT)*s*s);
  gc=1.0/ff+m.alt/XKMPER;
  gs=(1.0-FLAT)*(1.0-FLAT)/ff+m.alt/XKMPER;

  theta=gmst(mjd)+m.lng;
  dtheta=dgmst(mjd)*D2R/86400;

  pos->x=gc*cos(m.lat*D2R)*cos(theta*D2R)*XKMPER;
  pos->y=gc*cos(m.lat*D2R)*sin(theta*D2R)*XKMPER; 
  pos->z=gs*sin(m.lat*D2R)*XKMPER;
  vel->x=-gc*cos(m.lat*D2R)*sin(theta*D2R)*XKMPER*dtheta;
  vel->y=gc*cos(m.lat*D2R)*cos(theta*D2R)*XKMPER*dtheta; 
  vel->z=0.0;

  return;
}

// Solar position
void sunpos_xyz(double mjd,xyz_t *pos)
{
  double jd,t,l0,m,e,c,r;
  double n,s,ecl,ra,de;

  jd=mjd+2400000.5;
  t=(jd-2451545.0)/36525.0;
  l0=modulo(280.46646+t*(36000.76983+t*0.0003032),360.0)*D2R;
  m=modulo(357.52911+t*(35999.05029-t*0.0001537),360.0)*D2R;
  e=0.016708634+t*(-0.000042037-t*0.0000001267);
  c=(1.914602+t*(-0.004817-t*0.000014))*sin(m)*D2R;
  c+=(0.019993-0.000101*t)*sin(2.0*m)*D2R;
  c+=0.000289*sin(3.0*m)*D2R;

  r=1.000001018*(1.0-e*e)/(1.0+e*cos(m+c));
  n=modulo(125.04-1934.136*t,360.0)*D2R;
  s=l0+c+(-0.00569-0.00478*sin(n))*D2R;
  ecl=(23.43929111+(-46.8150*t-0.00059*t*t+0.001813*t*t*t)/3600.0+0.00256*cos(n))*D2R;

  ra=atan2(cos(ecl)*sin(s),cos(s));
  de=asin(sin(ecl)*sin(s));

  pos->x=r*cos(de)*cos(ra)*XKMPAU;
  pos->y=r*cos(de)*sin(ra)*XKMPAU;
  pos->z=r*sin(de)*XKMPAU;

  return;
}

// Compute precession angles
void precession_angles(double mjd0,double mjd,double *zeta,double *z,double *theta)
{
  double t0,t;

  // Time in centuries
  t0=(mjd0-51544.5)/36525.0;
  t=(mjd-mjd0)/36525.0;

  // Precession angles
  *zeta=(2306.2181+1.39656*t0-0.000139*t0*t0)*t;
  *zeta+=(0.30188-0.000344*t0)*t*t+0.017998*t*t*t;
  *zeta*=D2R/3600.0;
  *z=(2306.2181+1.39656*t0-0.000139*t0*t0)*t;
  *z+=(1.09468+0.000066*t0)*t*t+0.018203*t*t*t;
  *z*=D2R/3600.0;
  *theta=(2004.3109-0.85330*t0-0.000217*t0*t0)*t;
  *theta+=-(0.42665+0.000217*t0)*t*t-0.041833*t*t*t;
  *theta*=D2R/3600.0;
  
  return;
}

void compute_positions(char *tlefile,FILE *file,double ra0,double de0,double radius,int nmjd)
{
  int i,imode;
  orbit_t orb;
  FILE *fp=NULL;
  xyz_t satpos,satvel;
  double dx,dy,dz,r,ra,de,d,rsun,rearth;
  double psun,pearth,ptot;
  double a,b,c,age;
  char state[16];

  // Open TLE file
  fp=fopen(tlefile,"rb");
  if (fp==NULL) {
    printf("Error: TLE catalog %s not found, skipping\n",tlefile);
    return;
  }

  // Read TLEs
  while (read_twoline(fp,0,&orb)==0) {
    Isat=orb.satno;
    imode=init_sgdp4(&orb);
    
    // Skip on error
    if (imode==SGDP4_ERROR)
      continue;

    // Loop over times
    for (i=0;i<nmjd;i++) {
      // Satellite position
      satpos_xyz(p[i].mjd+2400000.5,&satpos,&satvel);

      // Check on radius
      r=sqrt(satpos.x*satpos.x+satpos.y*satpos.y+satpos.z*satpos.z);
      if (r>300000)
      	continue;
      
      // Relative to observer
      dx=satpos.x-p[i].obspos.x;  
      dy=satpos.y-p[i].obspos.y;
      dz=satpos.z-p[i].obspos.z;

      // Celestial position
      r=sqrt(dx*dx+dy*dy+dz*dz);
      ra=modulo(atan2(dy,dx),2.0*M_PI);
      de=asin(dz/r);
      
      // Correct for precession
      a=cos(de)*sin(ra+p[i].zeta);
      b=cos(p[i].theta)*cos(de)*cos(ra+p[i].zeta)-sin(p[i].theta)*sin(de);
      c=sin(p[i].theta)*cos(de)*cos(ra+p[i].zeta)+cos(p[i].theta)*sin(de);
      ra=modulo((atan2(a,b)+p[i].z)*R2D,360.0);
      de=asin(c)*R2D;
      
      // Check if nearby enough
      r=acos(sin(de0*D2R)*sin(de*D2R)+cos(de0*D2R)*cos(de*D2R)*cos((ra0-ra)*D2R))*R2D;
      if (r>radius)
	continue;

      // Satellite position relative to the Sun
      dx=-satpos.x+p[i].sunpos.x;  
      dy=-satpos.y+p[i].sunpos.y;
      dz=-satpos.z+p[i].sunpos.z;

      // Distances
      rsun=sqrt(dx*dx+dy*dy+dz*dz);
      rearth=sqrt(satpos.x*satpos.x+satpos.y*satpos.y+satpos.z*satpos.z);

      // Angles
      psun=asin(696.0e3/rsun)*R2D;
      pearth=asin(6378.135/rearth)*R2D;
      ptot=acos((-dx*satpos.x-dy*satpos.y-dz*satpos.z)/(rsun*rearth))*R2D;

      // Visibility state
      if (ptot-pearth<-psun) {
	strcpy(state,"eclipsed");
      } else if (ptot-pearth>-psun && ptot-pearth<psun) {
	strcpy(state,"umbra");
      } else if (ptot-pearth>psun) {
	strcpy(state,"sunlit");
      }

      // TLE age
      age=p[i].mjd+2400000.5-SGDP4_jd0;
      
      fprintf(file,"%05d,%014.8lf,%010.6lf,%+010.6lf,%s,%s,%.3f\n",orb.satno,p[i].mjd,ra,de,state,tlefile,age);
    }

  }
  fclose(fp);
  
  return;
}

// Present nfd
void nfd_now(char *s)
{
  time_t rawtime;
  struct tm *ptm;

  // Get UTC time
  time(&rawtime);
  ptm=gmtime(&rawtime);
    
  sprintf(s,"%04d-%02d-%02dT%02d:%02d:%02d",ptm->tm_year+1900,ptm->tm_mon+1,ptm->tm_mday,ptm->tm_hour,ptm->tm_min,ptm->tm_sec);
  
  return;
}

void usage()
{
  printf("satpredict t:l:c:R:D:r:L:B:H:o:\n\n");
  printf("t    date/time (yyyy-mm-ddThh:mm:ss.sss) [default: now]\n");
  printf("l    length (s) [default: 10s]\n");
  printf("n    number of points [default: 10]\n");
  printf("c    TLE catalog file [default: classfd.tle]\n");
  printf("R    R.A. (deg) [default: 0.0 deg]\n");
  printf("D    Decl. (deg) [default: 0.0 deg]\n");
  printf("r    radius (deg) [default: 10.0 deg]\n");
  printf("L    manual site longitude (deg) [default: 0.0 deg]\n");
  printf("B    manual site latitude (deg) [default: 0.0 deg]\n");
  printf("H    manual site elevation (m) [default: 0.0 deg]\n");
  printf("o    output csv file [default: results.csv]\n");
  printf("h    this help\n");

  return;
}

int main(int argc,char *argv[])
{
  int i,arg=0,nmjd=10,ntlefile=0;
  char nfd[32]="",outfile[LIM]="results.csv";
  char tlefile[NMAX][LIM];
  float length=10.0;
  double ra0=0.0,de0=0.0,radius=10.0;
  double mjd0,t;
  xyz_t obsvel;
  FILE *file;

  // Redirect stderr
  freopen("/tmp/stderr.txt","w",stderr);
  
  // Default options
  nfd_now(nfd);
  strcpy(tlefile[0],"classfd.tle");
  
  // Decode options
  if (argc>1) {
    while ((arg=getopt(argc,argv,"t:l:n:c:R:D:r:L:B:H:o:"))!=-1) {
      switch (arg) {
	
      case 't':
	strcpy(nfd,optarg);
	mjd0=nfd2mjd(nfd);
	break;
	
      case 'l':
	length=atof(optarg);
	break;

      case 'n':
	nmjd=atoi(optarg);
	break;
	
      case 'c':
	if (ntlefile==NMAX) {
	  printf("Error: Maximum number of TLE catalog files reached [%d]\n",NMAX);
	  return -1;
	}
	//	strcpy(tlefile[ntlefile++],optarg);
	snprintf(tlefile[ntlefile++],LIM-1,"%s",optarg);
	break;
	
      case 'o':
	strcpy(outfile,optarg);
	break;
	
      case 'R':
	ra0=(double) atof(optarg);
	break;
	
      case 'D':
	de0=(double) atof(optarg);
	break;
	
      case 'r':
	radius=(double) atof(optarg);
	break;
	
      case 'L':
	m.lng=atof(optarg);
	break;
	
      case 'B':
	m.lat=atof(optarg);
	break;
	
      case 'H':
	m.alt=atof(optarg)/1000.0;
	break;

      default:
	usage();
	return 0;
      }
    }
  } else {
    usage();
    return 0;
  }

  // Allocate
  p=(struct point *) malloc(sizeof(struct point)*nmjd);
  
  // Decode MJD  
  mjd0=nfd2mjd(nfd);

  // Initialize
  for (i=0;i<nmjd;i++) {
    // Compute time
    t=length*(float) i/(float) (nmjd-1);
    p[i].mjd=mjd0+t/86400.0;

    // Compute observer and sun position
    obspos_xyz(p[i].mjd,&p[i].obspos,&obsvel);
    sunpos_xyz(p[i].mjd,&p[i].sunpos);

    // Compute precession angles
    precession_angles(p[i].mjd,51544.5,&p[i].zeta,&p[i].z,&p[i].theta);
  }

  // Open output file
  file=fopen(outfile,"w");
  fprintf(file,"satno,mjd,ra,dec,state,tlefile\n");
  
  // Compute positions
  for (i=0;i<ntlefile;i++)
    compute_positions(tlefile[i],file,ra0,de0,radius,nmjd);

  // Close output file
  fclose(file);
  
  // Free
  free(p);

  // Close
  fclose(stderr);
  
  return 0;
}
