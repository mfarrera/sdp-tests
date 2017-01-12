
#include "vis.h"
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

// Visibility data
/*struct bl_data
{
    int antenna1, antenna2;
    int time_count;
    int freq_count;
    double *time;
    double *freq;
    double *uvw;
    double complex *vis;
    double complex *awkern;
};
struct vis_data
{
    int antenna_count;
    int bl_count;
    struct bl_data *bl;
};
*/

#define NUM_POL 1

double* arange(double start, double stop, double step){
	int len = (stop -start)/step;
	double *range = malloc(len*sizeof(double));
	double val=start;
	int i;
	for(i=0;i<len;i++){
		assert(val<=stop);
		range[i]=val;
		val=val+step;
	}
	return range;
}

double dot_product(double *a, double *b, int size){
 	double res=0.0;
	int i;
	for(i=0;i<size;i++){
		res+=a[i]*b[i];
	}
	return res;
}

void simulate_point(struct vis_data *vis, double l, double m){
    //Simulate visibilities for unit amplitude point source at
    //direction cosines (l,m) relative to the phase centre.

    //This includes phase tracking to the centre of the field (hence the minus 1
    //in the exponent.)

    //Note that point source is delta function, therefore the
    //FT relationship becomes an exponential, evaluated at
    //(uvw.lmn)

    //:param dist_uvw: :math:`(u,v,w)` distribution of projected baselines (in wavelengths)
    //:param l: horizontal direction cosine relative to phase tracking centre
    //:param m: orthogonal directon cosine relative to phase tracking centre

    // vector direction to source
    int i,j;
    double s[4];
    s[0] = l;
    s[1] = m;
    s[2] = sqrt(1.0-pow(l,2)-pow(m,2))-1.0;
    //s = numpy.array([l, m, numpy.sqrt(1 - l ** 2 - m ** 2) - 1.0])
    // complex valued Visibility data
    for(i=0;i< vis->bl_count;i++){
	for(j=0;j< vis->bl[i].time_count; j++){
		vis->bl[i].vis[j]+= cexp(-2i * M_PI * dot_product(&(vis->bl[i].uvw[j*3]),s,3));
    		}
	}
    //return numpy.exp(-2j * numpy.pi * numpy.dot(dist_uvw, s))
}



void xyz_to_uvw(struct antenna_configuration *ant, double ha, double dec){
    //Rotate :math:`(x,y,z)` positions in earth coordinates to
    //:math:`(u,v,w)` coordinates relative to astronomical source
    //position :math:`(ha, dec)`. Can be used for both antenna positions
    //as well as for baselines.

    //Hour angle and declination can be given as single values or arrays
    //of the same length. Angles can be given as radians or astropy
    //quantities with a valid conversion.

    //:param ant: :math:`(x,y,z)` co-ordinates of antennas in array
    // this is updated with the uvw coordinates for ha and dec
    //:param ha: hour angle of phase tracking centre (:math:`ha = ra - lst`)
    //:param dec: declination of phase tracking centre.
    
	
    int i;
    double u,v0,w,v;
    double x,y,z;
    
    for(i=0;i<ant->num_antennas; i++){
	x= ant->antenna[i].x;
	y= ant->antenna[i].y;
	z= ant->antenna[i].z;
    //# Two rotations:
    //#  1. by 'ha' along the z axis
    //#  2. by '90-dec' along the u axis
	u = x * cos(ha) - y * sin(ha);
	v0 = x* sin(ha) + y *cos(ha);
	w = z * sin(dec) -v0 *cos(dec);
	v = z * cos(dec) + v0 * sin(dec);
	ant->antenna[i].u=u;
	ant->antenna[i].v=v;
	ant->antenna[i].w=w;
    } 	

}




void baselines(struct antenna_configuration *ants_uvw,struct baselines_uvw *dist_uvw){

    //Compute baselines in uvw co-ordinate system from
    //uvw co-ordinate system station positions

    //:param ants_uvw: `(u,v,w)` co-ordinates of antennas in array
    // ants_uvw.shape[0] is num antenas len(ants_uvw)
    

    int i=0, j=0, k=0;
    int num_bl = ants_uvw->num_antennas * (ants_uvw->num_antennas-1)/2;
    dist_uvw->num_baselines = num_bl;
    dist_uvw->bl= malloc(num_bl*sizeof(struct uvw_bl));   
    for (i=0; i< ants_uvw->num_antennas; i++){
	for(j=i+1; j< ants_uvw->num_antennas; j++){
		dist_uvw->bl[k].antenna1=ants_uvw->antenna[i].antenna_id;
		dist_uvw->bl[k].antenna2=ants_uvw->antenna[j].antenna_id;
		dist_uvw->bl[k].u = ants_uvw->antenna[j].u - ants_uvw->antenna[i].u;
		dist_uvw->bl[k].v = ants_uvw->antenna[j].v - ants_uvw->antenna[i].v;
		dist_uvw->bl[k].w = ants_uvw->antenna[j].w - ants_uvw->antenna[i].w;
		k++;
	}
    }
	
}

void xyz_to_baseline(struct antenna_configuration *antennas, double *ha_range, double dec, struct vis_data *vis){
//dist_uvw = numpy.concatenate([baselines(xyz_to_uvw(ants_xyz, hax, dec)) for hax in ha_range])

	int i,j,k,t=0;
	double hax=0.0;
	struct baselines_uvw dist_uvw;
	vis->antenna_count=antennas->num_antennas;
	vis->bl_count = (vis->antenna_count)*(vis->antenna_count-1)/2;
	vis->bl = malloc(sizeof(struct bl_data)*vis->bl_count*LEN(ha_range)); 
	for(k=0;k<vis->bl_count;k++){
		vis->bl[k].time_count=LEN(ha_range);
		vis->bl[k].freq_count=1; //???
		vis->bl[k].time= malloc(vis->bl[k].time_count*sizeof(double));
		vis->bl[k].uvw= malloc(vis->bl[k].time_count*sizeof(double)*3);
		vis->bl[k].vis=malloc(vis->bl[k].time_count*vis->bl[k].freq_count*NUM_POL*sizeof(double complex));
		bzero(vis->bl[k].vis,vis->bl[k].time_count*vis->bl[k].freq_count*NUM_POL*sizeof(double complex));
	}
	for(i=0;i<LEN(ha_range);i++){
		hax=ha_range[i];
		xyz_to_uvw(antennas,hax,dec);
		baselines(antennas,&dist_uvw);	
		for (j=0;j<vis->bl_count;j++){
			vis->bl[j].time[i]=hax;
			vis->bl[j].uvw[i*3]= dist_uvw.bl[j].u;
			vis->bl[j].uvw[i*3+1]= dist_uvw.bl[j].v;
			vis->bl[j].uvw[i*3+2]= dist_uvw.bl[j].w;
		}
	}
	for(k=0;k<vis->bl_count;k++){
		vis->bl[k].antenna1=dist_uvw.bl[k].antenna1;
		vis->bl[k].antenna2=dist_uvw.bl[k].antenna2;
	}
	
	
}


//int generate_vis(char *vis_file, struct vis_data *vis,
//             double min_len, double max_len, double *ha_range, double dec){
int generate_vis(char *vis_file, struct vis_data *vis, double min_len, double max_len){

    //param vis_file has or we read it from array ants_xyz: :math:`(x,y,z)` co-ordinates of antennas in array
    //param ha_range: list of hour angle values for astronomical source as function of time
    //param dec: declination of astronomical source [constant, not :math:`f(t)`]
	struct antenna_configuration ac;
	if(vis_file ==NULL){
		// Read positions from the the array 
		ac.num_antennas = LEN(antenna_positions_VLA_a_hor)/3;
		ac.antenna = malloc(ac.num_antennas*sizeof(struct antenna_entry));
		int i,j=0;
		for(i=0;i<ac.num_antennas;i++){
			assert(j<LEN(antenna_positions_VLA_a_hor));
			ac.antenna[i].antenna_id=i;
			ac.antenna[i].x= antenna_positions_VLA_a_hor[j];
			j++;
			ac.antenna[i].y= antenna_positions_VLA_a_hor[j];
                        j++;
			ac.antenna[i].z= antenna_positions_VLA_a_hor[j];
                        j++;
		}
		xyz_to_baseline(&ac, arange(0,M_PI,0.04) , M_PI /4, vis); 
	}else{
		// Read the antenna configuration from a file
		// TODO
		fprintf(stderr,"Not implemented\n");
		exit(0);
	}
	simulate_point(vis,0.001,0.001);
	simulate_point(vis,0.0025,0.0025);
	return 0;

}

