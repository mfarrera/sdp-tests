

float Tobs=60*60*6;
float Tdump=0.1;

int Nf=64000; //Num frequency channels
int Ndump=toInt(ceil(Tobs/Tdump)); 
int Np =4; // Num polarizations 

int Ntime= 90; //Snapshot time
int Nb = 1; //Beams

float visibilityData[Nb][Nf][Ntime][Np][Ndump/Ntime];
float predictedVisibility[Nb][Nf][Ntime][Np][Ndump/Ntime];
float correctedVisibility[Nb][Nf][Ntime][Np][Ndump/Ntime];

typedef struct Component{
	int direction[3]; //type SkyCoord:
	float flux[Np][Nchannels];  
	float frequency[Nchannels];  //Number of channels? Nf??
	char shape[10]; // 'Point' or 'Gaussian'
	int params[10]; //shape dependent parameters
	char name[60]; 
}Component;	

Components locaSkyModel[Nb][MAXCOMP];

typedef struct voxel{
	unsigned char red,green,blue;
}Voxel;

Voxel skyImage[IMAGE_SIZE][IMAGE_SIZE];





