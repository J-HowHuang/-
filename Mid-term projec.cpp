//Mid-term Project
#include<iostream>
#include<cmath>
using namespace std;

double threatofP(double x0, double y0, int* x, int* y, double* r, int* p, int m);//threat of a point
double threatof(int** route, int* x, int* y, double* r, int* p, int m, int w, int t);//threat of a route
double length(int startX, int startY, int endX, int endY);
bool turnOrNot(int startX , int startY , int nowX , int nowY , int endX , int endY );
bool goThroughWall(int* x, int* y, double* r, int* p, int m ,int x0 , int y0) ;
bool lengthGoThroughWall(int* x, int* y, double* r, int* p, int m ,int x0, int y0, int x1, int y1) ;
void insertf(int **route,int x,int y, int endX,int endY) ;
void newRoute(int** route, int* x, int* y, double* r, int* p, int m, int w, int t);
const int MAX_CHANGING = 100;
int main(){
	int n = 0; //n: size of the map
	int m = 0; //m: number of threats
	int w = 0; //w: cost of changing direction
	int d = 0; //d: the limit of flying distance
	cin >> n >> m >> w >> d;
	
	int* x = new int [m]; //x[i]: the x-coordinate of the 'i'th threat
	int* y = new int [m]; //y[i]: the y-coordinate of the 'i'th threat
	double* r = new double [m]; //r[i]: the radius of the 'i'th threat
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
//algorithm start
    //remapping
    double pmax = 0;
    double *rp = new double [m];
    for(int i = 0 ; i < m ; i++ )
    {
        if(p[i] > pmax)
        {
            pmax = p[i];
        }
    }
    for(int i = 0 ; i < m ; i++ )
    {
        rp[i] =r[i] * (static_cast<double> (p[i]) / pmax ) * n / (n + m);
    }
    //open list
    int** open = new int* [n + 1];
    for(int i = 0; i < n + 1; i++){
        open[i] = new int [n + 1];
        for(int j = 0; j < n + 1; j++)
            open[i][j] = -1;
    }
    open[startX][startY] = 1;
    //source list
    int ***source = new int **[n + 1];
    double ***f = new double **[n + 1];
    for(int i = 0; i < n + 1; i++)
    {
        source[i] = new int *[n + 1];
        f[i] = new double *[n + 1];
        for(int j = 0; j < n + 1; j++){
            source[i][j] = new int[2];
            source[i][j][0] = startX;
            source[i][j][1] = startY;
            f[i][j] = new double [2];// 0:f,1:g
            f[i][j][0] = 9999999;
            f[i][j][1] = 0;
        }
    }
    f[startX][startY][0] = length(startX, startY, endX, endY);
    source[startX][startY][0] = -1;
    source[startX][startY][1] = -1;
    //
    while( open[endX][endY] != 0 )
    {
        int currentX = 0;
        int currentY = 0;
        double minf = INFINITY;
        for(int i = 0; i < n + 1; i++)//still need modified
            for(int j = 0; j < n + 1; j++)
                if(open[i][j] == 1)
                {
                cout << "(" << i << ", " << j << ") f: " << f[i][j][0] << " source: (" << source[i][j][0] << ", " << source[i][j][1] << ")\n";
                    if(f[i][j][0] < minf){
                        minf = f[i][j][0];
                        currentX = i;
                        currentY = j;
                    }
                }
      cout << "current position: (" << currentX << ", " <<currentY << ")\n";
        open[currentX][currentY] = 0;
        for(int i = -1; i <= 1; i += 1)
            for(int j = -1; j <= 1; j += 1)
            {
                //if it is not in neither open list or close list
                if(currentX + i < 0 || currentY + j < 0)
                    continue;
                if(currentX + i > n || currentY + j > n)
                    continue;
                if(open[currentX + i][currentY + j] != 0)
                {	
					if(abs(threatofP(currentX + i,currentY + j,x,y,rp,p,m)) > 0.0001 )
                    {
                        open[currentX + i][currentY + j] = 0;
                        continue;
                    }
                    //add it to open list
                    open[currentX + i][currentY + j] = 1;
                    source[currentX + i][currentY + j][0] = currentX;
                    source[currentX + i][currentY + j][1] = currentY;
                    
                    int nextsourceX = currentX + i;
                    int nextsourceY = currentY + j;
                    int finalsourceX = 0;
                    int finalsourceY = 0;
                    while( source[nextsourceX][nextsourceY][0] != -1 )
                    {
                        // checking wall or not
                        if(lengthGoThroughWall(x, y, rp, p, m, source[nextsourceX][nextsourceY][0], source[nextsourceX][nextsourceY][1], nextsourceX, nextsourceY))
                        {
                            int temp = nextsourceX;
                            nextsourceX = source[nextsourceX][nextsourceY][0];
                            nextsourceY = source[temp][nextsourceY][1];
                        }
                        else
                        {
                        	cout << "AA\n";
                            finalsourceX = source[nextsourceX][nextsourceY][0];
                            finalsourceY = source[nextsourceX][nextsourceY][1];
							nextsourceX = finalsourceX;
							nextsourceY = finalsourceY;
                        }
                    }
                    source[currentX + i][currentY + j][0] = finalsourceX;
                    source[currentX + i][currentY + j][1] = finalsourceY;
                    nextsourceX = currentX + i;
                    nextsourceY = currentY + j;
                    while(source[nextsourceX][nextsourceY][0] != -1){
                    	f[currentX + i][currentY + j][1] += length(nextsourceX, nextsourceY, source[nextsourceX][nextsourceY][0], source[nextsourceX][nextsourceY][1]);
                    	int tempX = source[nextsourceX][nextsourceY][0];
                    	int tempY = source[nextsourceX][nextsourceY][1];
						nextsourceX = tempX;
						nextsourceY = tempY;
					}
					f[currentX + i][currentY + j][0] = length(currentX + i, currentY + j, endX, endY) +  f[currentX + i][currentY + j][1];
                }
            }
    }
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
	for(int i = 1; i < turnCnt + 1; i++){
		cout << route[i + 1][0] << " " << route[i + 1][1] << " ";
	}
	cout << " risk: " << threatof(route, x, y, r, p, m, w, 0);
    

    return 0;
}
double threatofP(double x0, double y0, int* x, int* y, double* r, int* p, int m){
	double threat = 0;	
	for(int i = 0; i < m; i++)
		if((x0 - x[i]) * (x0 - x[i]) + (y0 - y[i]) * (y0 - y[i]) < r[i] * r[i])
			threat += p[i] * (r[i] - sqrt((x0 - x[i])*(x0 - x[i]) + (y0 - y[i])*(y0 - y[i])))/r[i];
	return threat;
}
double threatof(int** route, int* x, int* y, double* r, int* p, int m, int w, int t){
	//cout << "\nFunction is called.";//
	if(t != 0)
		newRoute(route, x, y, r, p, m, w, t);
	double leftLen = 0, threat = 0;
	int corner = 0;
	for(int i = 1; i <= (route[0][0] - t) || i == (route[0][0] + 1 - t); i++){
		//cout << "\n i is " << i;//
		if(i > 1){
			bool corTemp = turnOrNot(route[i - 1][0], route[i - 1][1], route[i][0], route[i][1], route[i + 1][0], route[i + 1][1]);
			corner += corTemp;
		}
		double len = length(route[i][0], route[i][1], route[i + 1][0], route[i + 1][1]);

		double cmpntX = (route[i + 1][0] - route[i][0]);//culculate x component 
		double cmpntY = (route[i + 1][1] - route[i][1]);
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
		if(i == (route[0][0] - t)){
			i = route[0][0];
			//cout << "\nroute[0][0] is " << route[0][0];//
		}
		//cout << "\ni is " << i;//
	}
	threat += w * corner;
	return threat;
}
void newRoute(int** route, int* x, int* y, double* r, int* p, int m, int w, int t){
	int k = route[0][0] - t;
	int** newRt = new int* [k];
	for(int i = 0; i <= k; i++){
		
	}
	
	// 
	threatof(newRt, x, y, r, p, m, w, t);
}
double length(int startX, int startY, int endX, int endY){
	double distance = 0;
	distance = (sqrt(pow((endX - startX), 2) + pow((endY - startY), 2)));
	return distance;
}
bool turnOrNot(int startX , int startY , int nowX , int nowY , int endX , int endY ){
	if(((nowY - startY) * (endX - nowX)) == ((endY - nowY) * (nowX - startX)) && (nowY - startY) * (endY - nowY) > 0)
		return true ;
	return false ;
}
bool goThroughWall(int* x , int* y , double* r , int* p , int m , int x0 , int y0){
	for(int i = 0 ; i < m ; i++)
	{
		double rTemp = pow((x[i] - x0),2) + pow((y[i] - y0),2);
		if(rTemp < pow(r[i] , 2))
		{
			return 1 ;
		}
	}
	return 0 ;
}
bool lengthGoThroughWall(int* x, int* y, double* r, int* p, int m ,int x0, int y0, int x1, int y1)
{
	for(int i = 0 ; i < m ; i++)
	{
		//A(x0 , y0),B(x1 , y1),P is the threat
		double cross = (x[i] - x0) * (x1 - x0) + (y[i] - y0) * (y1 - y0);//the inner product of AP and AB
		double distance = length(x0 , y0 , x1 , y1) ;//the distance between (x0 , y0) and (x1 , y1)
		double unit = cross / pow(distance , 2) ;
		double r2 = pow(r[i] , 2) ;
		//if the included angle of AB and AP is bigger the 90 degrees 
		if(cross <= 0) 
		{
			//the shortest distance between AB and P is the length of AP 
			double rTemp = pow((x[i] - x0),2) + pow((y[i] - y0),2) ;
			if(rTemp < r2){
				return 1 ;
			}
		}
		else if(cross >= length(x1, y1, x[i], y[i]))
		{
			//the shortest distance between AB and P is the length of AP 
			double rTemp = pow((x[i] - x1),2) + pow((y[i] - y1),2) ;
			if(rTemp < r2){
				return 1 ;
			}
		}
		//otherwise
		else{
			double xp = x0 + (x1 - x0) * unit ;//the x coordinate of the projection point of P on AB 
			double yp = y0 + (y1 - y0) * unit ;//the y coordinate of the projection point of P on AB
			double rTemp = pow((x[i] - xp),2) + pow((y[i] - yp),2) ; ;//the distance between the projection point and P 
			if(rTemp < r2){
				return 1 ;
			}
		}
	} 
	return 0 ;
}
void insertf(int **route,int x, int y, int endX, int endY){
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
