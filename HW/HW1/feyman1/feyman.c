#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>

# define A 3.0
# define B 2.0
# define C 1.0
# define DIM 3
# define H 0.001
# define STEPSZ sqrt((double)DIM * H)
# define DEFAULT_NK 6
# define DEFAULT_NI 16
# define DEFAULT_NJ 11

int i4_ceiling(double x)
{
  int value = (int)x;
  if (value < x)
    value = value + 1;
  return value;
}

int i4_min(int i1, int i2)
{
  int value;
  if (i1 < i2)
    value = i1;
  else
    value = i2;
  return value;
}

double potential(double a, double b, double c, double x, double y, double z)
{
  return 2.0 * (pow(x / 9.0, 2) + pow(y / 4.0, 2) + pow((double) z, 2)) + 1.25 + 1.0 / 9.0;
}

double r8_uniform_01(int *seed)
{
  int k;
  double r;

  k = *seed / 127773;

  *seed = 16807 * (*seed - k * 127773) - k * 2836;

  if (*seed < 0)
  {
    *seed = *seed + 2147483647;
  }
  r = (double)(*seed) * 4.656612875E-10;

  return r;
}

void timestamp(void)
{
#define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time(NULL);
  tm = localtime(&now);

  strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);

  printf("%s\n", time_buffer);

  return;
#undef TIME_SIZE
}

// print na stdout upotrebiti u validaciji paralelnog resenja
int main(int arc, char **argv)
{
  double chk;
  double dx;
  double dy;
  double dz;
  double err;
  int i;
  int j;
  int k;
  int n_inside;
  int seed = 123456789;
  int steps;
  int steps_ave;
  int trial;
  double us;
  double ut;
  double vh;
  double vs;
  double x;
  double x1;
  double x2;
  double x3;
  double y;
  double w;
  double w_exact;
  double we;
  double wt;
  double z;

  double time1, time2, elapsed_time;

  int N = atoi(argv[1]);
  timestamp();

  printf("A = %f\n", A);
  printf("B = %f\n", B);
  printf("C = %f\n", C);
  printf("N = %d\n", N);
  printf("H = %6.4f\n", H);

  err = 0.0;
  n_inside = 0;
  time1 = omp_get_wtime();
# pragma omp parallel for collapse(3) default(none) \
    shared(N) \
    reduction(+:err, n_inside) \
    private(i, j, k, x, y, z, w_exact, wt, steps, chk, trial, x1, x2, x3, w, ut, us, dx, dy, dz, vs, vh, we) \
    schedule(dynamic)
  for (i = 1; i <= DEFAULT_NI; i++)
  {
    for (j = 1; j <= DEFAULT_NJ; j++)
    {
      for (k = 1; k <= DEFAULT_NK; k++)
      {
        // x = ((double)(DEFAULT_NI - i) * (-A) + (double)(i - 1) * A) / (double)15;
        // y = ((double)(DEFAULT_NJ - j) * (-B) + (double)(j - 1) * B) / (double)10;
        // z = ((double)(DEFAULT_NK - k) * (-C) + (double)(k - 1) * C) / (double)5;
        x = (2.0 * i - 17.0) / 5.0;
        y = (2.0 * j - 12.0) / 5.0;
        z = (2.0 * k - 7.0) / 5.0;
        
        // Each thread gets a different seed, so that the random numbers are different.
        int local_seed = 123456789 + i * 10000 + j * 100 + k;

        chk = pow(x / A, 2) + pow(y / B, 2) + pow(z / C, 2);

        if (1.0 < chk)
        {
          // w_exact = 1.0;
          // wt = 1.0;
          // steps_ave = 0;
          // printf("  %7.4f  %7.4f  %7.4f  %10.4e  %10.4e  %10.4e  %8d\n",
          //        x, y, z, wt, w_exact, fabs(w_exact - wt), steps_ave);

          continue;
        }

        n_inside++;

        w_exact = exp(chk - 1.0);
        wt = 0.0;
        steps = 0;

        for (trial = 0; trial < N; trial++)
        {
          x1 = x;
          x2 = y;
          x3 = z;
          w = 1.0;
          chk = 0.0;
          while (chk < 1.0)
          {
            ut = r8_uniform_01(&local_seed);
            if (ut < 1.0 / 3.0)
            {
              us = r8_uniform_01(&local_seed) - 0.5;
              if (us < 0.0)
                dx = -STEPSZ;
              else
                dx = STEPSZ;
            }
            else
              dx = 0.0;

            ut = r8_uniform_01(&local_seed);
            if (ut < 1.0 / 3.0)
            {
              us = r8_uniform_01(&local_seed) - 0.5;
              if (us < 0.0)
                dy = -STEPSZ;
              else
                dy = STEPSZ;
            }
            else
              dy = 0.0;

            ut = r8_uniform_01(&local_seed);
            if (ut < 1.0 / 3.0)
            {
              us = r8_uniform_01(&local_seed) - 0.5;
              if (us < 0.0)
                dz = -STEPSZ;
              else
                dz = STEPSZ;
            }
            else
              dz = 0.0;

            vs = potential(A, B, C, x1, x2, x3);
            x1 = x1 + dx;
            x2 = x2 + dy;
            x3 = x3 + dz;

            steps++;

            vh = potential(A, B, C, x1, x2, x3);

            we = (1.0 - H * vs) * w;
            w = w - 0.5 * H * (vh * we + vs * w);

            chk = pow(x1 / A, 2) + pow(x2 / B, 2) + pow(x3 / C, 2);
          }
          wt += w;
        }
        wt /= (double)(N);
        // steps_ave = steps / (double)(N);

        err += pow(w_exact - wt, 2);

        // printf("  %7.4f  %7.4f  %7.4f  %10.4e  %10.4e  %10.4e  %8d\n",
        //        x, y, z, wt, w_exact, fabs(w_exact - wt), steps_ave);
      }
    }
  }
  time2 = omp_get_wtime();
  elapsed_time = time2 - time1;
  printf("\nElapsed time = %f seconds\n", elapsed_time);

  err = sqrt(err / (double)(n_inside));

  printf("\n\nRMS absolute error in solution = %e\n", err);
  timestamp();
  return 0;
}
