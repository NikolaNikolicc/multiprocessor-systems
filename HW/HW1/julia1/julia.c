# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>

# define DEFAULT_H 1000
# define DEFAULT_W 1000
# define DEFAULT_CNT 200
# define DEFAULT_FILENAME "julia"
# define DEFAULT_CI 0.156
# define DEFAULT_CR -0.8
# define DEFAULT_XL -1.5
# define DEFAULT_XR +1.5
# define DEFAULT_YB -1.5
# define DEFAULT_YT +1.5

int main (int argc, char *argv[]);
unsigned char *julia_set ( int w, int h, int cnt );
int julia ( float ar, float ai, int cnt );
void tga_write ( int w, int h, unsigned char rgb[], char *filename );
void timestamp ( );

int main (int argc, char *argv[] )  {
  int h = DEFAULT_H;
  int w = DEFAULT_W;
  int cnt = DEFAULT_CNT;
  char filename[256] = DEFAULT_FILENAME; 
  char buffer[256];
  unsigned char *rgb;
  double t1;
  double t2;

  if (argc == 4) {
  h = atoi(argv[1]);
  w = atoi(argv[2]);
  cnt = atoi(argv[3]);
  if (!h || !w || !cnt) return 1;
  }
  
  strcat(filename, "_");
  sprintf(buffer, "%d", h);
  strcat(filename, buffer);
  strcat(filename, "_");
  sprintf(buffer, "%d", w);
  strcat(filename, buffer);
  strcat(filename, "_");
  sprintf(buffer, "%d", cnt);
  strcat(filename, buffer);
  strcat(filename, ".tga");
    
  timestamp();
  printf ( "\n" );
  printf ( "JULIA Set\n" );
  printf ( "  Plot a version of the Julia set for Z(k+1)=Z(k)^2-0.8+0.156i\n" );

  rgb = julia_set ( w, h, cnt );

  tga_write ( w, h, rgb, filename );

  free ( rgb );

  printf ( "\n" );
  printf ( "JULIA set:\n" );
  printf ( "  Normal end of execution.\n" );

  timestamp();
  
  
  return 0;
}

unsigned char *julia_set ( int w, int h, int cnt )
{
  int i;
  int j;
  int juliaValue;
  int k;
  unsigned char *rgb;
  float x;
  float y;
  float factorX;
  float factorY;

  rgb = ( unsigned char * ) malloc ( w * h * 3 * sizeof ( unsigned char ) );

  for (i = 0; i < w * h * 3; i++) 
  {
    rgb[i] = 255;
  }

  k = 0;
  factorX = (float) (DEFAULT_XL - DEFAULT_XR) / (float) (w - 1);
  factorY = (float) (DEFAULT_YB - DEFAULT_YT) / (float) (h - 1);
  
  y = (float) DEFAULT_YB;
  
  for ( j = 0; j < h; j++ )
  { 
    x = (float) DEFAULT_XL;
    for ( i = 0; i < w; i++ )
    {
      juliaValue = julia ( x, y, cnt );
      if (juliaValue)
      {
        rgb[k]   = 0;
        rgb[k+1] = 0;
      }
      // k = (j * w + i) * 3;
      k += 3;
      // x = DEFAULT_XL - i * factorX;
      x -= factorX;
    }
    // y = DEFAULT_YB - j * factorY;
    y -= factorY;

  }
  return rgb;
}

int julia ( float ar, float ai, int cnt )
{
  int k;
  float t;

  for ( k = 0; k < cnt; k++ )
  {
    t  = ar * ar - ai * ai + DEFAULT_CR;
    ai = ar * ai + ai * ar + DEFAULT_CI;
    ar = t;

    if ( 1000 < ar * ar + ai * ai )
    {
      return 0;
    }
  }

  return 1;
}

void tga_write ( int w, int h, unsigned char rgb[], char *filename )
{
  FILE *file_unit;
  unsigned char header1[12] = { 0,0,2,0,0,0,0,0,0,0,0,0 };
  unsigned char header2[6] = { w%256, w/256, h%256, h/256, 24, 0 };

  file_unit = fopen ( filename, "wb" );

  fwrite ( header1, sizeof ( unsigned char ), 12, file_unit );
  fwrite ( header2, sizeof ( unsigned char ), 6, file_unit );

  fwrite ( rgb, sizeof ( unsigned char ), 3 * w * h, file_unit );

  fclose ( file_unit );

  printf ( "\n" );
  printf ( "TGA_WRITE:\n" );
  printf ( "  Graphics data saved as '%s'\n", filename );

  return;
}

void timestamp ( void )
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
