//Mid-term Project
#include<iostream>
#include<cmath>
using namespace std;

double threatofP(double x0, double y0, int* x, int* y, int* r, int* p, int m);//threat of a point
double threatof(int** route, int* x, int* y, int* r, int* p, int m, int w);//threat of a route
double length(int startX, int startY, int endX, int endY);
bool turnOrNot(int startX , int startY , int nowX , int nowY , int endX , int endY );
void insertf(int **route,int x,int y, int endX,int endY);
const int MAX_CHANGING = 10;
int main(){
	int n = 0; //n: size of the map
	int m = 0; //m: number of threats
	int w = 0; //w: cost of changing direction
	int d = 0; //d: the limit of flying distance
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
	cout << "input done\n";
//a* algorithm	
	//create openlist(1: open, 0: close, -1: not checked)
	int** open = new int* [n];
	for(int i = 0; i < n; i++){
		open[i] = new int [n];
		for(int j = 0; j < n; j++)
			open[i][j] = -1;
	}
	
	cout << "openlist done\n";
//add start point to open list	
	int openCnt = 1;
	open[startX][startY] = 1;
	
	//f[][] is the min. approx cost of a point
	double** f = new double* [n];
	for(int i = 0; i < n; i++){
		f[i] = new double [n];
		for(int j = 0; j < n; j++)
			f[i][j] = INFINITY;
	}
	int** originRoute = new int* [3];
	for(int i = 0; i < 3; i++){
		originRoute[i] = new int [2];
	}
	originRoute[0][0] = 0;
	originRoute[0][1] = 0;
	originRoute[1][0] = startX;
	originRoute[1][1] = startY;
	originRoute[2][0] = endX;
	originRoute[2][1] = endY;
	f[startX][startY] = threatof(originRoute, x, y, r, p, m, w);
	cout << threatof(originRoute, x, y, r, p, m, w) << "\n";
	
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
	cout << "process 1 done\n";
	while(openCnt != 0){
		//search the point with min. approx cost in the open list
		int currentX;
		int currentY;
		double minf = 99999999;
		for(int i = 0; i < n; i++)//still need modified
			for(int j = 0; j < n; j++)
				if(open[i][j] == 1)
					if(f[i][j] < minf){
						minf = f[i][j];
						currentX = i;
						currentY = j;
					}
 		cout << "find the min in open\n";
		//add it into close list
		open[currentX][currentY] = 0;
		openCnt -= 1;
		//for each point near it
		for(int i = -1; i <= 1; i++)
			for(int j = -1; j <= 1; j++){
				if(currentX + i < 0 || currentY + j < 0)
					continue;
				else if(currentX + i > n || currentY + j > n)
					continue;
				//if it is not in neither open list or close list
				if(open[currentX + i][currentY + j] == -1){
					//add it to open list
					open[currentX + i][currentY + j] = 1;
					openCnt += 1;
					source[currentX + i][currentY + j][0] = currentX;
					source[currentX + i][currentY + j][1] = currentY;
					//calculate the f(t)
					int nextsourceX = currentX + i;
					int nextsourceY = currentY + j;
					int turnCnt = 0;
					int** tempRoute = new int* [MAX_CHANGING + 3];
					for(int k = 0; k < MAX_CHANGING + 3; k++){
						tempRoute[k] = new int[2];
						tempRoute[k][0] = 0;
						tempRoute[k][1] = 0;
					}
					cout << "A";
					//trace back from t and save the route in reverse	
					tempRoute[1][0] = endX;
					tempRoute[1][1] = endY;
					while(nextsourceX != -1){
						tempRoute[turnCnt + 2][0] = nextsourceX;
						tempRoute[turnCnt + 2][1] = nextsourceY;
						turnCnt += 1;
						tempRoute[0][0] += 1;
						tempRoute[0][1] += 1;
						int tempX = source[nextsourceX][nextsourceY][0];
						int tempY = source[nextsourceX][nextsourceY][1];
						nextsourceX = tempX;
						nextsourceY = tempY;
					}
					cout << "B";
					turnCnt -= 1;
					tempRoute[0][0] -= 1;
					tempRoute[0][1] -= 1;
					
					//reverse the temp route to get the route in correct order
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
					cout << "C";
					//if point t has better performance pass through former turn point straightly, update the source
					for(int k = 0; k < turnCnt; k++){
						cout << k;
						cout << threatof(route, x, y, r, p, m, k);
						if(threatof(route, x, y, r, p, m, k) < f[currentX + i][currentY + j]){
							f[currentX + i][currentY + j] = threatof(route, x, y, r, p, m, k);
							source[currentX + i][currentY + j][0] = route[turnCnt - k][0];
							source[currentX + i][currentY + j][1] = route[turnCnt - k][1];
						}
						cout << k;
					}
					cout << 3*i+j+5 << "\n";
					
				}
				//if it is in open list 
				else if(open[currentX + i][currentY + j] == 1){
					//consider?
				}
			}
		if(open[endX][endY] == 0)
			break;
	}
	//trace back from t and save the route in reverse
	int turnCnt = 0;	
	int nextsourceX = source[endX][endY][0];
	int nextsourceY = source[endX][endY][1];
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
		int tempX = source[nextsourceX][nextsourceY][0];
		int tempY = source[nextsourceX][nextsourceY][1];
		nextsourceX = tempX;
		nextsourceY = tempY;
	}
	turnCnt -= 1;
	tempRoute[0][0] -= 1;
	tempRoute[0][1] -= 1;
		
	//reverse the temp route to get the route in correct order
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
	
	cout << turnCnt << " ";
	for(int i = 0; i < turnCnt + 1; i++){
		cout << route[i + 1][0] << " " << route[i + 1][1] << " ";
	}
	cout << route[turnCnt + 2][0] << " " << route[turnCnt + 2][1] << " ";
	cout << threatof(route, x, y, r, p, m, 0);
	return 0;
	
}
double threatofP(double x0, double y0, int* x, int* y, int* r, int* p, int m){
	double threat = 0;
	
	for(int i = 0; i < m; i++)
		if((x0 - x[i]) * (x0 - x[i]) + (y0 - y[i]) * (y0 - y[i]) < r[i] * r[i])
			threat += p[i] * (r[i] - sqrt((x0 - x[i])*(x0 - x[i]) + (y0 - y[i])*(y0 - y[i])))/r[i];
}
double threatof(int** route, int* x, int* y, int* r, int* p, int m, int w){
	double leftLen = 0, threat = 0;
	int corner = 0;
	for(int i = 1; i <= (route[0][0] + 1); i++){
		if(i > 1){
			bool corTemp = turnOrNot(route[i - 1][0], route[i - 1][1], route[i][0], route[i][1], route[i + 1][0], route[i + 1][1]);
			corner += corTemp;
		}
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
			threat += threatofP(tempX, tempY, x, y, r, p, m);
		}
		leftLen = 0;
	}
	threat += w * corner;
	return threat;
}
double length(int startX, int startY, int endX, int endY){
	double distance = 0;
	distance = (sqrt(pow((endX - startX), 2) + pow((endY - startY), 2)));
	return distance;
}
bool turnOrNot(int startX , int startY , int nowX , int nowY , int endX , int endY ){
	if((nowY - startY / nowX - startX  ) == (endY - nowY / endX - nowX ))
		return true ;
	return false ;
}
void insertf(int **route,int x,int y, int endX,int endY){
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
