//Mid-term Project
#include<iostream>
#include<cmath>
using namespace std;

double threatof(double x0, double y0, int* x, int* y, int* r, int* p, int m);//threat of a point
double threatof(int** route, int* x, int* y, int* r, int* p, int m);//threat of a route
double approxCost(int x0, int y0, int* x, int* y, int* r, int* p, int m);
int main(){
	int n = 0; //n: size of the map
	int m = 0; //m: number of threats
	int w = 0; //w: cost of changing direction
	int d = 0; //d: the limit of flying distance
	const int MAX_CHANGING = 10;
	cin >> n >> m >> w >> d;
	
	int* x = new int [m]; //x[i]: the x-coordinate of the 'i'th threat
	int* y = new int [m]; //y[i]: the y-coordinate of the 'i'th threat
	int* r = new int [m]; //r[i]: the radius of the 'i'th threat
	int* p = new int [m]; //p[i]: the power of the 'i'th threat
	for(int i = 0; i < m; i++)
		cin >> x[i];
	for(int i = 0; i < m; i++)
		cin >> y[i];
	for(int i = 0; i < m; i++)
		cin >> r[i];
	for(int i = 0; i < m; i++)
		cin >> p[i];

	int startX = 0, startY = 0; //the start point (startX, startY)
	int endX = 0, endY = 0; //the end point (endX, endY)
	cin >> startX >> startY >> endX >> endY;
//a* algorithm	
	//create openlist(1: open, 0: close, -1: not checked)
	int** open = new int* [n];
	for(int i = 0; i < n; i++){
		open[i] = new int [n];
		for(int j = 0; j < n; j++)
			open[i][j] = -1;
	}
	
	
//add start point to open list	
	int openCnt = 1;
	open[startX][startY] = 1;
	
	//f[][] is the min. approx cost of a point
	double** f = new double* [n];
	for(int i = 0; i < n; i++){
		f[i] = new double [n];
		for(int j = 0; j < n; j++)
			f[i][j] = -1;
	}
	f[startX][startY] = threatof(route, x, y, r, p);
	
	//source[][][0] is the x coordinate of the source of the point, while 1 is y
	int*** source = new int** [n];
	for(int i = 0; i < n; i++){
		source[i] = new int*[n];
		for(int j = 0; j < n; j++){
			source[i][j] = new int[2];
			source[i][j][0] = startX;
			source[i][j][1] = startY;
		}
	}
	source[startX][startY][0] = -1;
	source[startX][startY][1] = -1;
	
	while(openCnt != 0){
		//search the point with min. approx cost in the open list
		int currentX;
		int currentY;
		double minf = 9999999
		for(int i = 0; i < n; i++)//still need modified
			for(int j = 0; j < n; j++)
				if(open[i][j] == 1)
					if(f[i][j] < minf && f[i][j] > 0){
						minf = f[i][j];
						currentX = i;
						currentY = j;
					}
 
		//add it into close list
		open[currentX][currentY] = 0;
		//for each point near it
		for(int i = -1, i <= 1; i++)
			for(int j = -1; j <= 1; j++){
				//if it is not in neither open list or close list
				if(open[currentX + i][currentY + j] == -1){
					//add it to open list
					open[currentX + i][currentY + j] = 1;
					openCnt += 1;
					source[currentX + i][currentY + j][0] = currentX;
					source[currentX + i][currentY + j][1] = currentY;
					//calculate the f(t)
					double minf = 9999999;
					int nextsourceX = currentX + i;
					int nextsourceY = currentY + j;
					int turnCnt = 0;
					int** tempRoute = new int* [MAX_CHANGING + 3];
					for(int k = 0; k < MAX_CHANGING + 3; k++){
						tempRoute[k] = new int[2];
						tempRoute[k][0] = 0;
						tempRoute[k][1] = 0;
					}
						
					tempRoute[1][0] = endX;
					tempRoute[1][1] = endY;
					while(nextsourceX != -1){
						tempRoute[turnCnt + 2][0] = nextsourceX;
						tempRoute[turnCnt + 2][1] = nextsourceY;
						turnCnt += 1;
						tempRoute[0][0] += 1;
						tempRoute[0][1] += 1;
						tempX = source[nextsourceX][nextsourceY][0];
						tempY = source[nextsourceX][nextsourceY][1];
						nextsourceX = tempX;
						nextsourceY = tempY;
					}
					turnCnt -= 1;
					tempRoute[0][0] -= 1;
					tempRoute[0][1] -= 1;
					int** route = new int* [MAX_CHANGING + 3];
					for(int k = 0; k < MAX_CHANGING + 3; k++){
						route[k] = new int [2];
						route[k][0] = 0;
						route[k][1] = 0;
					}
					route[0][0] = tempRoute[0][0];
					route[0][1] = tempRoute[0][1];
					for(int k = 0; k < turnCnt + 2; k++){
						route[k + 1][0] = tempRoute[turnCnt + 2 - k][0];
						route[k + 1][1] = tempRoute[turnCnt + 2 - k][1];
					}
					for(int k = 0; 
				}
				//if it is in open list 
				else if(open[currentX + i][currentY + j] == 1){
					//consider?
				}
			}
		if(open[endX][endY] == 0)
			break;
	}
	
	return 0;
	
}
double threatof(double x0, double y0, int* x, int* y, int* r, int* p, int m){
	double threat = 0;
	
	for(int i = 0; i < m; i++)
		if((x0 - x[i]) * (x0 - x[i]) + (y0 - y[i]) * (y0 - y[i]) < r[i] * r[i])
			threat += p[i] * (r[i] - sqrt((x0 - x[i])*(x0 - x[i]) + (y0 - y[i])*(y0 - y[i])))/r[i];
}
double threatof(int** route, int* x, int* y, int* r, int* p){
	
}
double approxCost(int x0, int y0, int* x, int* y, int* r, int* p, int m){
    
	
}
void insertf(int **route,int x,int y, int endX,int enddY)
{
    for(int i = 0 ; i < MAX_CHANGING + 3 ; i++ )
    {
        if( ( route[i][0] == endX ) && ( route[i][1] == endY) )
        {
            route[i + 1][0] = route[i][0];
            route[i + 1][1] = route[i][1];
            route[i][0] = x;
            route[i][1] = y;
            break;
        }
    }
    
}
