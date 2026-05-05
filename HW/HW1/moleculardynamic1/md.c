# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include <math.h>

int main ( int argc, char *argv[] );
void compute ( int np, int nd, double pos[], double vel[], 
  double f[], double *pot, double *kin );
double cpu_time ( );
double dist ( int nd, double r1[], double r2[], double dr[] );
void initialize ( int np, int nd, double pos[], double vel[], double acc[] );
void r8mat_uniform_ab ( int m, int n, int *seed, double r[] );
void timestamp ( );
void update ( int np, int nd, double pos[], double vel[], double f[], 
  double acc[], double dt );

int main ( int argc, char *argv[] )
{
  double *acc;
  double ctime;
  double dt;
  double e0;
  double *force;
  int id;
  double kinetic;
  int nd;
  int np;
  double *pos;
  double potential;
  int step;
  int step_num;
  int step_print;
  double *vel;

  timestamp ( );
  printf ( "\n" );
  printf ( "MD\n" );
  printf ( "  C version\n" );
  printf ( "  A molecular dynamics program.\n" );
/*
  Get the spatial dimension.
*/
  if ( 1 < argc )
  {
    nd = atoi ( argv[1] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter ND, the spatial dimension (2 or 3).\n" );
    scanf ( "%d", &nd );
  }
//
//  Get the number of particles.
//
  if ( 2 < argc )
  {
    np = atoi ( argv[2] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter NP, the number of particles (500, for instance).\n" );
    scanf ( "%d", &np );
  }
//
//  Get the number of time steps.
//
  if ( 3 < argc )
  {
    step_num = atoi ( argv[3] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter STEP_NUM, the number of time steps (500 or 1000, for instance).\n" );
    scanf ( "%d", &step_num );
  }
//
//  Get the time steps.
//
  if ( 4 < argc )
  {
    dt = atof ( argv[4] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter DT, the size of the time step (0.1, for instance).\n" );
    scanf ( "%lf", &dt );
  }
/*
  Report.
*/
  printf ( "\n" );
  printf ( "  ND, the spatial dimension, is %d\n", nd );
  printf ( "  NP, the number of particles in the simulation, is %d\n", np );
  printf ( "  STEP_NUM, the number of time steps, is %d\n", step_num );
  printf ( "  DT, the size of each time step, is %f\n", dt );
/*
  Allocate memory.
*/
  acc = ( double * ) malloc ( nd * np * sizeof ( double ) );
  force = ( double * ) malloc ( nd * np * sizeof ( double ) );
  pos = ( double * ) malloc ( nd * np * sizeof ( double ) );
  vel = ( double * ) malloc ( nd * np * sizeof ( double ) );
/*
  This is the main time stepping loop:
    Compute forces and energies,
    Update positions, velocities, accelerations.
*/
  printf ( "\n" );
  printf ( "  At each step, we report the potential and kinetic energies.\n" );
  printf ( "  The sum of these energies should be a constant.\n" );
  printf ( "  As an accuracy check, we also print the relative error\n" );
  printf ( "  in the total energy.\n" );
  printf ( "\n" );
  printf ( "      Step      Potential       Kinetic        (P+K-E0)/E0\n" );
  printf ( "                Energy P        Energy K       Relative Energy Error\n" );
  printf ( "\n" );

  step_print = 0;
  
  ctime = cpu_time ( );

  initialize ( np, nd, pos, vel, acc );
  compute ( np, nd, pos, vel, force, &potential, &kinetic );
  e0 = potential + kinetic;
  printf ( "  %8d  %14f  %14f  %14e\n", step, potential, kinetic,
      ( potential + kinetic - e0 ) / e0 );
    step_print += step_num / 10;
  for ( step = 1; step <= step_num; step++ )
  {
    update ( np, nd, pos, vel, force, acc, dt );

    compute ( np, nd, pos, vel, force, &potential, &kinetic );

    if ( step == step_print )
    {
      printf ( "  %8d  %14f  %14f  %14e\n", step, potential, kinetic,
       ( potential + kinetic - e0 ) / e0 );
      step_print += step_num / 10;
    }

  }
/*
  Report timing.
*/
  ctime = cpu_time ( ) - ctime;
  printf ( "\n" );
  printf ( "  Elapsed cpu time: %f seconds.\n", ctime );
/*
  Free memory.
*/
  free ( acc );
  free ( force );
  free ( pos );
  free ( vel );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "MD\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}

void compute ( int np, int nd, double pos[], double vel[], 
  double f[], double *pot, double *kin )
{
  double d;
  double d2;
  int j;
  int k;
  double ke;
  double pe;
  double PI2 = 3.141592653589793 / 2.0;
  double rij[3];

  pe = 0.0;
  ke = 0.0;

  for ( k = 0; k < np; k++ )
  {
/*
  Compute the potential energy and forces.
*/
    f[k*nd] = 0.0;
    f[1+k*nd] = 0.0;
    if (nd == 3) f[2+k*nd] = 0.0;

    for ( j = 0; j < np; j++ )
    {
      if ( k != j )
      {
        d = dist ( nd, pos+k*nd, pos+j*nd, rij );
/*
  Attribute half of the potential energy to particle J.
*/
        if ( d < PI2 )
        {
          d2 = d;
        }
        else
        {
          d2 = PI2;
        }

        pe = pe + 0.5 * pow ( sin ( d2 ), 2 );

        f[k*nd] = f[k*nd] - rij[0] * sin ( 2.0 * d2 ) / d;
        f[1+k*nd] = f[1+k*nd] - rij[1] * sin ( 2.0 * d2 ) / d;
        if (nd == 3) f[2+k*nd] = f[2+k*nd] - rij[2] * sin ( 2.0 * d2 ) / d;
      }
    }
/*
  Compute the kinetic energy.
*/
    ke = ke + vel[k*nd] * vel[k*nd];
    ke = ke + vel[1+k*nd] * vel[1+k*nd];
    if (nd == 3) ke = ke + vel[2+k*nd] * vel[2+k*nd];
  }

  ke = ke * 0.5;
  
  *pot = pe;
  *kin = ke;

  return;
}

double cpu_time ( )
{
  double value;

  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

  return value;
}
/******************************************************************************/

double dist ( int nd, double r1[], double r2[], double dr[] )
{
  double d;
  double tmp;
  // 0
  d = 0.0;
  tmp = r1[0] - r2[0];
  dr[0] = tmp;
  d += tmp * tmp;
  // 1
  tmp = r1[1] - r2[1];
  dr[1] = tmp;
  d += tmp * tmp;
  if (nd == 3)
  {
    // 2
    tmp = r1[2] - r2[2];
    dr[2] = tmp;
    d += tmp * tmp;
  }
  d = sqrt ( d );

  return d;
}
/******************************************************************************/

void initialize ( int np, int nd, double pos[], double vel[], double acc[] )
{
  int j;
  int seed;
/*
  Set positions.
*/
  seed = 123456789;
  r8mat_uniform_ab ( nd, np, &seed, pos );
/*
  Set velocities.
*/
  for ( j = 0; j < np; j++ )
  {
    vel[j*nd] = 0.0;
    vel[1+j*nd] = 0.0;
    if (nd == 3) vel[2+j*nd] = 0.0;
  }
/*
  Set accelerations.
*/
  for ( j = 0; j < np; j++ )
  {
    acc[j*nd] = 0.0;
    acc[1+j*nd] = 0.0;
    if (nd == 3) acc[2+j*nd] = 0.0;
  }

  return;
}

void r8mat_uniform_ab ( int m, int n, int *seed, double r[] )
{
  int i;
  int j;
  int k;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_UNIFORM_AB - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + 2147483647;
      }
      r[i+j*m] = 46.56612875E-10 * ( double ) ( *seed );
    }
  }

  return;
}

void timestamp ( )
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}

void update ( int np, int nd, double pos[], double vel[], double f[], 
  double acc[], double dt )
{
  int j;
  int position;

  for ( j = 0; j < np; j++ )
  {
    position = j * nd;
    pos[position] = pos[position] + vel[position] * dt + 0.5 * acc[position] * dt * dt;
    vel[position] = vel[position] + 0.5 * dt * ( f[position] + acc[position] );
    acc[position] = f[position];

    position = 1 + j * nd;
    pos[position] = pos[position] + vel[position] * dt + 0.5 * acc[position] * dt * dt;
    vel[position] = vel[position] + 0.5 * dt * ( f[position] + acc[position] );
    acc[position] = f[position];

    if (nd == 3)
    {
      position = 2 + j * nd;
      pos[position] = pos[position] + vel[position] * dt + 0.5 * acc[position] * dt * dt;
      vel[position] = vel[position] + 0.5 * dt * ( f[position] + acc[position] );
      acc[position] = f[position];
    }
  }

  return;
}