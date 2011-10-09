/*Wiener filter*/

#include <stdio.h>
#include "imlib-1.9.15/Imlib/Imlib.h"
//#include <complex.h>
#include <fftw3.h>
//#include <rfftw.h>
#include <math.h>
#include <stdlib.h>

#include "wiener.h"

#define TINY 1.0e-20

fftw_complex *out;

double *
genpsf (double *psf, int w, int h, ImlibImage * impsf, ImlibData * idpsf,
    double sigma)
{
  int i, j;
  double sum = 0.0;
  double max = 0.0;
  psf = malloc (sizeof (double) * w * h);
  for (j = 0; j < h; j++)
    {
      for (i = 0; i < w; i++)
    {
      psf[i + w * j] =
        (1.0 / (sqrt (2.0 * M_PI) * sigma))
            * exp (-0.5 * (pow (i - w / 2, 2) + pow (j - h / 2, 2)) / pow (sigma, 2));
      sum += psf[i + w * j];
      if (psf[i + w * j] > max)
        {
          max = psf[i + w * j];
        }
    }
    }
  printf ("Som is %f \n", sum);
  printf ("Max is %f \n", max);
  for (j = 0; j < h; j++)
    {
      for (i = 0; i < w; i++)
    {
      int index = (i + w * j) * 3;
      psf[i + w * j] = psf[i + w * j] / sum;
    }
    }
  sum = 0;
  max = 0;
  for (j = 0; j < h; j++)
    {
      for (i = 0; i < w; i++)
    {
      sum += psf[i + w * j];
      if (psf[i + w * j] > max)
        {
          max = psf[i + w * j];
        }
    }
    }
  printf ("Som is %f \n", sum);

  return psf;

}

/*extract from an RGB image the R, G of B channel*/
double *
kanaal (int w, int h, ImlibImage * im, int chan)
{
  int i, j;
  int index;
  double *out = malloc (sizeof (double) * w * h);

  for (i = 0; i < w; i++)
    {
      for (j = 0; j < h; j++)
    {
      index = (i + w * j) * 3;
      out[i + w * j] = im->rgb_data[index + chan];
    }
    }
  return out;
}

double *
regroup (int w, int h, fftw_complex * convred, double *flimoutn)
{
/*regroup the four quadrants after the DFT*/
  int i, j;

  for (j = 0; j < h; j++)
    {
      for (i = 0; i < w; i ++)
    {
      int k, l;
      int index1, index2;
      double temp;

      index1 = i + w * j;
      if (i < w / 2)
        {
          k = i + (w / 2);
        }
      if (j < h / 2)
        {
          l = j + (h / 2);
        }
      if (i >= w / 2)
        {
          k = i - (w / 2);
        }
      if (j >= h / 2)
        {
          l = j - (h / 2);
        }
      index2 = k + w * l;
      flimoutn[index2] = sqrt (((convred[index1][0] * convred[index1][0]
                      +  convred[index1][1] * convred[index1][1])));

    }
    }
  return flimoutn;
}

/*norm PSF to 1*/
double *
normpsf (double *psfch, int w, int h)
{
  int i, j;
  double sum;
  sum = 0.0;
  for (j = 0; j < h; j++)
    {
      for (i = 0; i < w; i++)
    {
      sum += psfch[i + w * j];
    }
    }
  for (j = 0; j < h; j++)
    {
      for (i = 0; i < w; i++)
    {
      psfch[i + w * j] = psfch[i + w * j] / sum;
    }
    }
  sum = 0.0;
  for (j = 0; j < h; j++)
    {
      for (i = 0; i < w; i++)
    {
      sum += psfch[i + w * j];
    }
    }
  return psfch;
}

/*convolve two images*/
fftw_complex *
wiener (double *kanaal, double *kanpsf, int w, int h, fftw_complex * out,
      double K)
{
  fftw_complex *im, *fim;
  fftw_complex *psfin, *fpsf, *conv;
  fftw_plan p, pinv;
  int i, j;
  /*input image*/
  im    = (fftw_complex *) malloc (w * h * sizeof (fftw_complex));
  /*fourier transf. of input image*/
  fim   = (fftw_complex *) malloc (w * h * sizeof (fftw_complex));
  psfin = (fftw_complex *) malloc (w * h * sizeof (fftw_complex));
  fpsf  = (fftw_complex *) malloc (w * h * sizeof (fftw_complex));
  conv  = (fftw_complex *) malloc (w * h * sizeof (fftw_complex));
  for (j = 0; j < h; j++)
    {
      for (i = 0; i < w; i++)
    {
      im[i + j * w][0] = kanaal[i + j * w];
      im[i + j * w][1] = 0;
    }
    }
  for (j = 0; j < h; j++)
    {
      for (i = 0; i < w; i++)
    {
      psfin[i + j * w][0] = kanpsf[i + j * w];
      psfin[i + j * w][1] = 0;
    }
    }
  p = fftw_plan_dft_2d (h, w, &im[0], &fim[0], FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  p = fftw_plan_dft_2d (h, w, &psfin[0], &fpsf[0], FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);

//   fftwnd_one (p, &im[0], &fim[0]);
//   fftwnd_one (p, &psfin[0], &fpsf[0]);
  for (j = 0; j < h; j++)
    {
      for (i = 0; i < w; i++)
    {
      conv[i + j * w][0] =
        (fim[i + j * w][0] * fpsf[i + j * w][0] +
         fim[i + j * w][1] * fpsf[i + j * w][1])
         / (w * h * (K + pow (fpsf[i + j * w][0], 2) +
                (pow (fpsf[i + j * w][1], 2))));
      conv[i + j * w][1] =
        (fim[i + j * w][0] * fpsf[i + j * w][1] -
         fim[i + j * w][1] * fpsf[i + j * w][0])
         / (w * h * (K + pow(fpsf[i + j * w][0], 2) +
                (pow(fpsf[i + j * w][1], 2))));
    }
    }
  pinv = fftw_plan_dft_2d (h, w, &conv[0], &out[0], FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(pinv);
  //fftwnd_one (pinv, &conv[0], &out[0]);
  free (im);
  free (fim);
  free (psfin);
  free (fpsf);
  free (conv);
  return out;
}


int wFilter (const char* inFile, const char* outFile, double sigma, double K)
{

  int teller, maxiter, chan;
  double *in, *est, *psf, *a, *b, *estblur;
  fftw_complex *cestblur, *cest;

  int w, h, i, j, iter;
  double total;

  ImlibData  *id, *idpsf, *idout;
  ImlibImage *im, *impsf, *imout;

  Display *disp;

  disp = XOpenDisplay (NULL);
  id = Imlib_init (disp);
  idpsf = Imlib_init (disp);
  idout = Imlib_init (disp);

  im    = Imlib_load_image (id, inFile);
  impsf = Imlib_load_image (idpsf, inFile);
  imout = Imlib_load_image (idout, inFile);

  /* Suck the image's original width and height out of the Image structure */
  w = im->rgb_width;
  h = im->rgb_height;
  out = (fftw_complex *) malloc (w * h * sizeof (fftw_complex));
  a = malloc (sizeof (double) * w * h);
  b = malloc (sizeof (double) * w * h);
  est = malloc (sizeof (double) * w * h);
  estblur = malloc (sizeof (double) * w * h);

  /*in */
  in = kanaal (w, h, im, 0);

  total = 0;
  for (j = 0; j < h; j++)
    {
      for (i = 0; i < w; i++)
    {
      total += in[i + w * j];
    }
    }
  printf (" \t \tvoor g \t %f \n", total);

  /*generate gaussian psf with given sigma*/
  genpsf (psf, w, h, impsf, idpsf, sigma);

  /*norm PSF to 1 */
  normpsf (psf, w, h);


  total = 0.0;
  for (j = 0; j < h; j++)
    {
      for (i = 0; i < w; i++)
    {
      total += psf[i + w * j];
    }
    }
  printf (" \t \tvoor psf \t %f \n", total);

  /*wienerfilter image given the PSF */
  wiener (in, psf, w, h, cestblur, K);
  regroup (w, h, cestblur, estblur);

  /*check: is total intensity of the filtered image the same? */
  total = 0.0;
  for (j = 0; j < h; j++)
    {
      for (i = 0; i < w; i++)
    {
      total += estblur[i + w * j];
    }
    }
  printf ("\t \t voor (f*H) \t %f \n", total);


/* OUTPUT */

  for (j = 0; j < h; j++)
    {
      for (i = 0; i < w; i++)
    {
      int index = (i + j * w) * 3;
      int index2 = i + w * j;
      int index3 = (w-1-i + (h-1-j) * w) * 3;
      if (estblur[index2] < 0)
        {
          estblur[index2] = 0;
        }
      if (estblur[index2] > 255)
        {
          estblur[index2] = 255;
        }
      imout->rgb_data[index3] = (int) estblur[index2];
      imout->rgb_data[index3 + 1] = (int) estblur[index2];
      imout->rgb_data[index3 + 2] = (int) estblur[index2];
    }
    }
  Imlib_save_image (idout, imout, outFile, NULL);

  return 1;
}
 
