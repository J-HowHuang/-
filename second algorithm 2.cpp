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
const int MAX_CHANGING = 1000;
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
    // 演算法死大頭！
    //全能地圖改造王
    int pmax = 0;
    int *rp = new int [m];
    for(int i = 0 ; i < m ; i++ )
    {
        if(p[i] > pmax)
        {
            pmax = p[i];
            
        }
    }
    for(int i = 0 ; i < m ; i++ )
    {
        rp[i] =r[i] * ( p[i] / pmax );
    }
    //open力死頭！
    int** open = new int* [n + 1];
    for(int i = 0; i < n + 1; i++){
        open[i] = new int [n + 1];
        for(int j = 0; j < n + 1; j++)
            open[i][j] = -1;
    }
    //source
    int ***source = new int **[n + 1];
    int ***f = new int **[n + 1];
    for(int i = 0; i < n + 1; i++)
    {
        source[i] = new int *[n + 1];
        f[i] = new int *[n + 1];
        for(int j = 0; j < n + 1; j++){
            source[i][j] = new int[2];
            source[i][j][0] = startX;
            source[i][j][1] = startY;
            f[i][j] = new int [2];// 0:f,1:g
            f[i][j][0] = 0;
            f[i][j][1] = 0;
        }
    }
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
                    //    cout << "(" << i << ", " << j << ") f: " << f[i][j] << " source: (" << source[i][j][0] << ", " << source[i][j][1] << ")\n";
                    if(f[i][j][0] < minf){
                        minf = f[i][j][0];
                        currentX = i;
                        currentY = j;
                    }
                }
        open[currentX][currentY] = 0;
        for(int i = -unit; i <= unit; i += unit)
            for(int j = -unit; j <= unit; j += unit)
            {
                //if it is not in neither open list or close list
                //    cout << "(" << currentX + i << ", " <<currentY + j << ")\n";
                if(currentX + i < 0 || currentY + j < 0)
                    continue;
                if(currentX + i > n || currentY + j > n)
                    continue;
                if(open[currentX + i][currentY + j] != 0)
                {
                    //add it to open list
                    open[currentX + i][currentY + j] = 1;
                    source[currentX + i][currentY + j][0] = currentX;
                    source[currentX + i][currentY + j][1] = currentY;
                    if(threatofP(currentX,currentY,x,y,p,m) != 0 )
                    {
                        open[currentX][currentY] = 0;
                        continue;
                    }
                    int nextsourceX = currentX;
                    int nextsourceY = currentY;
                    int finalsourceX = 0;
                    int finalsourceY = 0;
                    while( nextsourceX != -1 )
                    {
                        // checking wall or not
                        if(blocked)
                        {
                            int temp = nextsourceX;
                            nextsourceX = source[nextsourceX][nextsourceY][0];
                            nextsourceY = source[temp][nextsourceY][1];
                        }
                        else
                        {
                            int temp = nextsourceX;
                            finalsourceX = nextsourceX;
                            finalsourceY = nextsourceY;
                            nextsourceX = source[nextsourceX][nextsourceY][0];
                            nextsourceY = source[temp][nextsourceY][1];
                            
                        }
                    }
                    source[currentX + i][currentY + j][0] = finalsourceX;
                    source[currentX + i][currentY + j][1] = finalsourceY;
                    
                    
                }
            }
        
    }
    

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
        if( i > ( route[0][0] + 1 - t )  )
        {
            continue;
        }
        //
        
        if( ( i == route[0][0] + 1 - t ) )
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
