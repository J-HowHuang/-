//Mid-term Project
#include<iostream>
#include<cmath>
#include<ctime>
using namespace std;

double threatofP(double x0, double y0, int* x, int* y, int* r, int* p, int m);//threat of a point
double threatof(int** route, int* x, int* y, int* r, int* p, int m, int w, int t);//threat of a route
double length(int startX, int startY, int endX, int endY);
bool turnOrNot(int startX , int startY , int nowX , int nowY , int endX , int endY );
void insertf(int **route,int x,int y, int endX,int endY);
const int MAX_CHANGING = 10000;
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
//a* algorithm	
	clock_t s = clock();
	//set the unit
	int unit = 0;
	if(n < 100)
		unit = 1;
	else if(n < 400)
		unit = 3;
        else if(n < 750)
                unit = 8;
	else
		unit = 97 ;
	unit = 21;
	startX = startX + (endX - startX) % unit;
	startY = startY + (endY - startY) % unit;
	//create openlist(1: open, 0: close, -1: not checked)
	int** open = new int* [n + 1];
	for(int i = 0; i < n + 1; i++){
		open[i] = new int [n + 1];
		for(int j = 0; j < n + 1; j++)
			open[i][j] = -1;
	}
	
//add start point to open list	
	int openCnt = 1;
	open[startX][startY] = 1;
	
	//f[][] is the min. approx cost of a point
	double** f = new double* [n + 1];
	for(int i = 0; i < n + 1; i++){
		f[i] = new double [n + 1];
		for(int j = 0; j < n + 1; j++)
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

	f[startX][startY] = threatof(originRoute, x, y, r, p, m, w, 0);
//	cout << threatof(originRoute, x, y, r, p, m, w, 0) << "\n";

	
	//source[][][0] is the x coordinate of the source of the point, while 1 is y
	int*** source = new int** [n + 1];
	for(int i = 0; i < n + 1; i++){
		source[i] = new int*[n + 1];
		for(int j = 0; j < n + 1; j++){
			source[i][j] = new int[2];
			source[i][j][0] = startX;
			source[i][j][1] = startY;
		}
	}
	source[startX][startY][0] = -1;
	source[startX][startY][1] = -1;
	while(true){
		//search the point with min. approx cost in the open list
		int currentX = 0;
		int currentY = 0;
		double minf = INFINITY;
		for(int i = 0; i < n + 1; i++)//still need modified
			for(int j = 0; j < n + 1; j++)
				if(open[i][j] == 1){
				//	cout << "(" << i << ", " << j << ") f: " << f[i][j] << " source: (" << source[i][j][0] << ", " << source[i][j][1] << ")\n";
					if(f[i][j] < minf){
						minf = f[i][j];
						currentX = i;
						currentY = j;
					}
				}
					
 	//	cout << "current position: (" << currentX << ", " << currentY << ")\n";
		//add it into close list
		open[currentX][currentY] = 0;
		openCnt -= 1;
		//for each point near it
		for(int i = -unit; i <= unit; i += unit)
			for(int j = -unit; j <= unit; j += unit){
				//if it is not in neither open list or close list
				if(currentX + i < 0 || currentY + j < 0)
					continue;
				if(currentX + i > n || currentY + j > n)
					continue;
				if(open[currentX + i][currentY + j] != 0){
					//add it to open list
					open[currentX + i][currentY + j] = 1;
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
					//if point t has better performance pass through former turn point straightly, update the source
					for(int k = 0; k < turnCnt; k++){
						double threatofRoute = threatof(route, x, y, r, p, m, w, k);
					//	cout << "if ignore " << k << " corner " << threatofRoute << "\n";
						if(threatofRoute < f[currentX + i][currentY + j]){
							double lng = (length(currentX + i, currentY + j, endX, endY));
							f[currentX + i][currentY + j] = threatofRoute;
							source[currentX + i][currentY + j][0] = route[turnCnt - k][0];
							source[currentX + i][currentY + j][1] = route[turnCnt - k][1];
						}
					}
					
					for(int k = 0; k < MAX_CHANGING + 3; k++){
						delete[] tempRoute[k];
						delete[] route[k];
					}
					delete[] tempRoute;
					delete[] route;
				}
				//if it is in open list 
				
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
    //refine
    int correctturnCnt = 0;
    bool cor = 0;
    bool* route1 = new bool [MAX_CHANGING + 3];
    for(int k = 0; k < MAX_CHANGING + 3; k++)
    {
        route1[k] = 0;
    }
    for(int i = 2 ; i < turnCnt + 2 ; i++ )
    {
        cor = turnOrNot(route[i-1][0], route[i-1][1], route[i][0], route[i][0], route[i+1][0], route[i+1][1]);
        if(cor)
        {
            correctturnCnt++;
            route1[i] = 1;
        }
    }
    
    
    //
    cout << correctturnCnt << " ";
    for(int i = 2; i < turnCnt + 2; i++)
    {
        if(route1[i] == 0)
        {
            continue;
        }
        cout << route[i][0] << " " << route[i][1] << " ";
    }
	/*
	cout << turnCnt << " ";
	for(int i = 1; i < turnCnt + 1; i++){
		cout << route[i + 1][0] << " " << route[i + 1][1] << " ";
	}
     */
//	cout << "risk: " << threatof(route, x, y, r, p, m, w, 0);
	clock_t e = clock();
//	cout << "\ntime = " << e - s;
	return 0;
	
}
double threatofP(double x0, double y0, int* x, int* y, int* r, int* p, int m){
	double threat = 0;
	
	for(int i = 0; i < m; i++)
		if((x0 - x[i]) * (x0 - x[i]) + (y0 - y[i]) * (y0 - y[i]) < r[i] * r[i])
			threat += p[i] * (r[i] - sqrt((x0 - x[i])*(x0 - x[i]) + (y0 - y[i])*(y0 - y[i])))/r[i];
	return threat;
}
double threatof(int** route, int* x, int* y, int* r, int* p, int m, int w, int t){
    double leftLen = 0, threat = 0;
    int corner = 0;
    double len = 0;
    double cmpntX = 0;//culculate x component
    double cmpntY = 0;
    for(int i = 1; i <= ( route[0][0] + 1 ) ; i++){
        if( (i > ( route[0][0] - t ) ) && ( i < ( route[0][0] + 1 ) ) )
        {
            continue;
        }
        //
        
        if( ( i == route[0][0] - t ) )
        {
            if(i > 1)
            {
                bool corTemp = turnOrNot(route[i - 1][0], route[i - 1][1], route[i][0], route[i][1], route[route[0][0] + 1][0], route[route[0][0] + 1][1]);
                corner += corTemp;
            }
            len = length(route[i][0], route[i][1], route[route[0][0] + 1][0], route[route[0][0] + 1][1]);
            cmpntX = (route[route[0][0]+ 1][0] - route[i][0]);//culculate x component
            cmpntY = (route[route[0][0] + 1][1] - route[i][1]);
        }
        else if (i == route[0][0] + 1)
        {
            if(i > 1)
            {
                bool corTemp = turnOrNot(route[route[0][0] - t][0], route[route[0][0] - t][1], route[i][0], route[i][1], route[i + 1][0], route[i + 1][1]);
                corner += corTemp;
            }
            len = length(route[i][0], route[i][1], route[i + 1][0], route[i + 1][1]);
            cmpntX = (route[i + 1][0] - route[i][0]);//culculate x component
            cmpntY = (route[i + 1][1] - route[i][1]);
        }
        else
        {
            if(i > 1)
            {
                bool corTemp = turnOrNot(route[i - 1][0], route[i - 1][1], route[i][0], route[i][1], route[i + 1][0], route[i + 1][1]);
                corner += corTemp;
            }
            len = length(route[i][0], route[i][1], route[i + 1][0], route[i + 1][1]);
            cmpntX = (route[i + 1][0] - route[i][0]);//culculate x component
            cmpntY = (route[i + 1][1] - route[i][1]);
        }
        
        //
        
        
        double tempX = route[i][0];
        double tempY = route[i][1];
        if(abs(leftLen) > 0.0001){
            tempX += cmpntX / len * (1 - leftLen);
            tempY += cmpntY / len * (1 - leftLen);
            threat += threatofP(tempX, tempY, x, y, r, p, m);
        }
        int intLen = static_cast<int>(len - leftLen);
        leftLen = len - intLen - leftLen;
        for(int k = 0; k < intLen ; k++){
            tempX += cmpntX / len;
            tempY += cmpntY / len;
            threat += threatofP(tempX, tempY, x, y, r, p, m);
        }
        //cout << "\ni is " << i;//
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
	if(((nowY - startY) * (endX - nowX)) == ((endY - nowY) * (nowX - startX)) && (nowY - startY) * (endY - nowY) >= 0)
		return false;
	return true;
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
