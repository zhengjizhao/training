#ifndef GRAPH1010
#define GRAPH1010
#include <stdio.h>
#include <stdlib.h>
#include <vector>

using namespace std;
#define GRAPH_ENABLED true

class plot {
public:
 FILE *gp;
 bool enabled,persist;
 plot(bool _persist=false,bool _enabled=GRAPH_ENABLED) {
 enabled=_enabled;
 persist=_persist;
  if (enabled) {
	if(persist)
		  gp=popen("gnuplot -persist","w");
	else
		  gp=popen("gnuplot","w");
  }
 }

 void plot_data(vector<double> x,const char* style="points",const char* title="Data") {
 if(!enabled)
	return;
  fprintf(gp,"set title '%s' \n",title);
  fprintf(gp,"plot '-' w %s \n",style);
  for(int k=0;k<x.size();k++) {
   fprintf(gp,"%f\n",x[k]);
  }
  fprintf(gp,"e\n");
  fflush(gp);
 }

void plot_data(vector<double>& x,vector<double>& y,double max_x,double max_y, double min_x,double min_y,double time, const char* style="points") {
    if(!enabled)
	    return;
	fprintf(gp,"set title '%.20E (seconds)' \n",time);
    fprintf(gp,"set xrange [%.20E : %.20E] \n",min_x,max_x);
    fprintf(gp,"set yrange [%.20E : %.20E] \n",min_y,max_y);
    fprintf(gp,"set format y '%s' \n","%.10E");
    fprintf(gp,"plot '-' w %s \n",style);
    int N_points = 1000;
    for(int k=0;k<x.size();k+=x.size()/N_points) {
        fprintf(gp,"%.20E %.20E \n",x[k],y[k]);
    }
    fprintf(gp,"e\n");
    fflush(gp);
   }
   ~plot() {
   if(enabled)
       pclose(gp);
   }
 
};

/*
int main(int argc,char **argv) {
 plot p;
 for(int a=0;a<100;a++) {
 vector<double> x,y;
 for(int k=a;k<a+200;k++) {
   x.push_back(k);
   y.push_back(k*k);
 }
  p.plot_data(x,y);
 }
 return 0;
}
*/

#endif
