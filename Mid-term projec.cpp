//Mid-term Project
#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

/*
// PSEUDOCODE OF ThreatOfRoute
	
	
*/

double threatof(double x0, double y0, int* x, int* y, int* r, int* p, int m);//threat of a point
double threatof(int** route, int* x, int* y, int* r, int* p, int m);//threat of a route
double length(int startX, int startY, int endX, int endY);
double approxCost(int x0, int y0, int* x, int* y, int* r, int* p, int m);
const int MAX_N = 1000;

int main(){
	int n = 0; // n: size of the map
	int m = 0; // m: number of threats
	int w = 0; // w: cost of changing direction
	int d = 0; // d: the limit of flying 
	const int MAX_CHANGING = 10;
	cin >> n >> m >> w >> d;
	
	int* x = new int [m]; // x[i]: the x-coordinate of the 'i'th threat
	int* y = new int [m]; // y[i]: the y-coordinate of the 'i'th threat
	int* r = new int [m]; // r[i]: the radius of the 'i'th threat
	int* p = new int [m]; // p[i]: the power of the 'i'th threat
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
	
	//approxCost[i][j] saves the approxCost of point (i, j)
	double** approxCost = new double* [n];
	for(int i = 0; i < n; i++)
		approxCost[i] = new double [n];
	
	//route have yet selected
	int route[MAX_CHANGING + 3][2] = {0};
	route[1][0] = startX;
	route[1][1] = startY;
	route[2][0] = endX;
	route[2][1] = endY;
	// start point is the 1th point.
	// route[i][0] indicates the x coordinary of the ith point
	// route[i][1] indicates the y coordinary of the ith point
	// route[0][0] indicates the number of curve points (start point and terminal is not included) 
	
	int openCnt = 1;
	open[startX][startY] = 1;
	while(openCnt != 0){
		//search the point with min. approx cost in the open list
		int currentX;
		int currentY;
		//add it into close list
		open[currentX][currentY] = 0;
		//for each point near it
		for(int i = -1; i <= 1; i++)
			for(int j = -1; j <= 1; j++){
				//if it is not in neither open list or close list
				if(open[currentX + i][currentY + j] == -1){
					//add it to open list
					open[currentX + i][currentY + j] = 1;
					//record the approx cost?
				}
				//if it is in open list
				else if(open[currentX + i][currentY + j] == 1){
					//consider?
				}
			}
	}
	
	return 0;
	
}
double threatof(double x0, double y0, int* x, int* y, int* r, int* p, int m){
	double threat = 0;
	
	for(int i = 0; i < m; i++)
		if((x0 - x[i]) * (x0 - x[i]) + (y0 - y[i]) * (y0 - y[i]) < r[i] * r[i])
			threat += p[i] * (r[i] - sqrt((x0 - x[i])*(x0 - x[i]) + (y0 - y[i])*(y0 - y[i])))/r[i];
}
double threatof(int** route, int* x, int* y, int* r, int* p, int m){
	double leftLen = 0, threat = 0;
	for(int i = 1; i <= (route[0][0] + 1); i++){
		double len = length(route[i][0], route[i][1], route[i + 1][0], route[i + 1][1]);
		double cmpntX = sqrt(pow((route[i + 1][0] - route[i][0]), 2));
		double cmpntY = sqrt(pow((route[i + 1][1] - route[i][1]), 2));
		int intLen = static_cast<int>(len);
		leftLen = len - intLen;
		double tempX = route[i][0];
		double tempY = route[i][1];
		for(int k = 1; k < intLen ; k++){
			tempX += cmpntX / len;
			tempY += cmpntY / len;
			threat += threatof(tempX, tempY, x, y, r, p, m);
		}
		leftLen = 0;
	}
	return threat;
}
double length(int startX, int startY, int endX, int endY){
	double distance = 0;
	distance = (sqrt(pow((endX - startX), 2) + pow((endY - startY), 2)));
	return distance;
}
double approxCost(int x0, int y0, int* x, int* y, int* r, int* p, int m){
	
}
