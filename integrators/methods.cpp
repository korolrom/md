#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum eMethod{
  eEuler,
  eVerlet,

  eNum,
};

int main (int argc, const char * argv[]) {
  eMethod method   = eEuler;
  int     step_num = 60;
  if (argc < 2){ // parse args
    printf("usage:integration verlet 60\n");
    exit(-1);
  }else{
    if (strncmp(argv[1], "euler" , 6) == 0) { method = eEuler;  }
    if (strncmp(argv[1], "verlet", 6) == 0) { method = eVerlet; }
    sscanf(argv[2], "%d", &step_num);
  }
  double dt = 1.0 / (double(step_num) / 2.0);
  double t  = 0.0;
  double x  = 0.0;
  double a  = 0.0;
  double v  = 0.0;
  double m  = 10.0;
  double k  = 1000.0;
  double sx = 1.0;
  double f  = 0.0f; // F = k(x-sx)

  printf("# time position\n");
  //  printf("%f %f\n", t, x);
  switch(method) {
  case eEuler:
    for(int i = 0; i < step_num; i++) {
      f = -k * (x - sx);
      a = f / m;
      x += v * dt;
      v += a * dt;
      printf("%f %f\n", t, x);
      t += dt;
    }       
    break;
  case eVerlet:
    double prev_x = 0.0;
    for(int i = 0; i < step_num; i++) {
      f = -k * (x - sx);
      a = f / m;
      double tmp = x;
      x += (x - prev_x) + a * dt * dt;
      prev_x = tmp;
      printf("%f %f\n", t, x);
      t += dt;
     }
     break;
  }
  return 0;
}
