/* Final Assignment - Matekane Airstrip
**created December 8, 2011 by metafalsevac
**
*/
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

#define PI 3.14159265//defines PI for use in converting to radians  (not actually necessary as gamma is always in radians)
#define T0 518.67 //sea level temperature in Rankine 
#define P0 2116.224 //sea level pressure in lb/ft^2
#define a -0.003566 //Temperature Lapse Rate in Rankine/ft
#define R 1716.0 //Specific Gas Constant for air in ft*lbf/slug/degree Rankine
#define W 2712.5 //modified weight in pounds (-12.5% of original weight)
#define m 84.2391304348 //modified mass in slugs calculated from W/g=3100/32.2 (reduced by 12.5%)
#define S /*174//*/195.75 //modified reference area in sq ft (+12.5%)
#define CLMax 2.195 //CLMax
#define T0SL /*635.5//*/714.9375 //modified max sea level thrust in pounds (+12.5%)
#define aSL 0.0334 //airspeed variation coefficient for max thrust at sea level in lbs*sec^2/ft^2
#define g 32.2 //gravitational acceleration in ft/sec^2
#define Vxw0 83.0 //initial condition for Vxw
#define X0 0.0 //initial x distance
#define gamma0 0.0 //intitial flight path angle gamma
#define H0 7649.63 //initial altitude in feet

using namespace std;

double Ground (double x); //calculates elevation of ground at horizontal distance x
double Press (double h); //calculates pressure at altitude
double Temp (double h); //calculates temperature at altitude
double Rho (double h); //calculates rho at altitude
double Veas (double vtrue, double h); //calculates equivalent airspeed
double Vdot (double vxw, double gamma, double h, double CD); //Vxwind derivative
double Gammadot (double vxw, double gamma, double h, double CL); //gamma derivative
double Hdot (double vxw, double gamma); //H derivative
double Xdot (double vxw, double gamma); //X derivative
void RK (double vxw, double gamma, double h, double x, double y[6], double CD, double CL); //function for 4th order runge kutta integration
                                                                                           //this function is void because it needs to return 5 values
                                                                                           //it does so using arrays, which are treated as pointers
int main(){
double y[6], vxw, veas, gamma, CLmin, CLmax, step, h, x;//declare array used to return results of integration and variables for I.C.s
double CL = 2.195;
double CD = 0.025 + 0.056*CL*CL;

ofstream outfile ("output.txt");//open output file

//initial conditions
vxw = 83.0;
veas = 0.0; //this will be updated once RK begins
gamma = 0.0;
h = 7649.63;
x = 0.0;

cout << "Takeoff From Matekane Airstrip\n\n";
cout << "This program will simulate takeoff from the Matekane Airstrip \nat various lift coefficients.\n\nInput minimum lift coefficient: ";
cin >> CLmin;
cout << "\nInput maximum lift coefficient: ";
cin >> CLmax;
cout << "\nInput lift coefficient step size: ";
cin >> step;


//populate y array with initial conditions. this is really just for show since it is done again within for loop
y[0] = vxw; //velocity in the x-wind direction
y[1] = gamma; //flight path angle
y[2] = h; //altitude of plane
y[3] = x; //horizontal x-earth location of plane 
y[4] = 7649.63; //ground elevation
y[5] = 0.0; //equivalent airspeed

//run RK for multiple CL values with this for loop
for (CL = CLmin; CL <= CLmax; CL += step) {
y[0] = vxw;//must reset initial conditions each iteration
y[1] = gamma;
y[2] = h;
y[3] = x;
y[4] = 7649.63;
CD = 0.025 + 0.056*CL*CL;

    while (y[3]<=10000.0) {//while still in the valley, the while loop actually runs the integration past the far end of the valley to determine if the plane continues to gain altitude
      
      //outfile << "Vxw = " << y[0] << " gamma = " << y[1] << " h = " << y[2] << " x = " << y[3] << " ground = " << y[4] << endl;
      outfile << y[3] << "     " << y[2] << "     " << y[0] << "      " << y[4] << "     " << y[1] << "     " << y[5] << "     " << CL << endl;
      
      //this is the important line of this while loop. the output lines are for debugging or creating plots
      RK (y[0], y[1], y[2], y[3], y, CD, CL);//call runge kutta function
      
      
      //cout << "2 Vxw = " << y[0] << " gamma = " << y[1] << " h = " << y[2] << " x = " << y[3] << " ground = " << y[4] << endl << endl;
      //outfile << "2 Vxw = " << y[0] << " gamma = " << y[1] << " h = " << y[2] << " x = " << y[3] << " ground = " << y[4] << endl << endl;
      
      /*//optimize Cl         //this optimization code is commented out because it results in the airplane doing loops, which is obviously undesirable 
      if (95.0 > y[5]){
               CL -= .01;
               }
      if (y[5] > 107.0) {
               CL += .01;
               }*/
      
      
      if ((y[4]>=y[2])){//did plane crash?
         cout << "CL = " << CL << " CRASH!!!" << endl;
         cout << "Vxw = " << y[0] << " gamma = " << y[1] << " Alt = " << y[2] << " x = " << y[3] << " Elev = " << y[4] << "\nEquiv Airspeed = " << y[5] << endl << endl;//coordinates and flight conditions at crash
         break;
                      }//if
         
                           }//while
      if ((y[3]>6233.58)&&y[2]>y[4]){
         cout << "CL = " << CL << " SURVIVE!!!" << endl;
         cout << "Vxw = " << y[0] << " gamma = " << y[1] << " Alt = " << y[2] << " x = " << y[3] << " Elev = " << y[4] << "\nEquiv Airspeed = " << y[5] << endl << endl;
                     }//if
                                       }//for


outfile.close();//close output file    
cin.get();
cin.get();
return 0;    
}

double Ground (double x) {//determines ground elevation at horizontal x distance
       double H = -7.54673e-12*x*x*x*x + 9.16044e-8*x*x*x - 1.8723e-4*x*x - 5.91546e-1*x + 7.64963e3;
       return H;
       }

double Press (double h) {//determine pressure at altitude 
       double p = P0 * pow((1 + (a*h/T0)), 5.2561);
       return p;
       }
       
double Temp (double h) {//determine temp at altitude
       double t = T0 + a*h;
       return t;
       }
       
double Rho (double h) {//determine dens at altitude
       double rho = Press(h)/(R*Temp(h));
       return rho;
       }
       
double Veas (double vtrue, double h) {
       double veas = pow((Rho(h)/Rho(0)),0.5)*vtrue;
       return veas;
       }
       
double Vdot (double vxw, double gamma, double h, double CD) {//derivative of velocity i.e. acceleration
       double vdot = ( (Rho(h)/(m*Rho(0))) * (T0SL-(aSL*vxw*vxw)) - (Rho(h)*vxw*vxw*S*CD/(2*m)) - g*sin(gamma) );
       return vdot;
       }
       
double Gammadot (double vxw, double gamma, double h, double CL) {//derivative of flight path angle
       double gammadot = -( (g*cos(gamma)/vxw) - (Rho(h)*vxw*S*CL/(2*m)) );
       return gammadot;
       }
       
double Hdot (double vxw, double gamma) {//verticle speed
       double hdot = vxw*sin(gamma);
       return hdot;
       }
       
double Xdot (double vxw, double gamma) {//horizontal speed
       double xdot = vxw*cos(gamma);
       return xdot;
       }
       
void RK (double vxw, double gamma, double h, double x, double y[], double CD, double CL) {//runge kutta. this function is void and uses Arrays to return results of integration
       double delta = 0.1;//increment
       double k11, k12, k13, k14, k21, k22, k23, k24, k31, k32, k33, k34, k41, k42, k43, k44;//first index corresponds to each k of RK4, 2nd index corresponds to each y equation (i.e. Vdot,
       //k1 values for each of the 4 y functions                                            //Gammadot, Hdot, or Xdot)
       k11 = Vdot (vxw, gamma, h, CD);//k1 of Vdot
       k12 = Gammadot (vxw, gamma, h, CL);//k1 of Gammadot
       k13 = Hdot (vxw, gamma);//k1 of Hdot
       k14 = Xdot (vxw, gamma);//k1 of Xdot
       //k2 values
       k21 = Vdot ((vxw+.5*k11*delta), (gamma+.5*k12*delta), (h+.5*k13*delta), CD);//k2 of Vdot
       k22 = Gammadot ((vxw+.5*k11*delta), (gamma+.5*k12*delta), (h+.5*k13*delta), CL);//k2 of Gammadot
       k23 = Hdot ((vxw+.5*k11*delta), (gamma+.5*k12*delta));//k2 of Hdot
       k24 = Xdot ((vxw+.5*k11*delta), (gamma+.5*k12*delta));//k2 of Xdot
       //k3 values
       k31 = Vdot ((vxw+.5*k21*delta), (gamma+.5*k22*delta), (h+.5*k23*delta), CD);//and so on and so forth...
       k32 = Gammadot ((vxw+.5*k21*delta), (gamma+.5*k22*delta), (h+.5*k23*delta), CL);
       k33 = Hdot ((vxw+.5*k21*delta), (gamma+.5*k22*delta));
       k34 = Xdot ((vxw+.5*k21*delta), (gamma+.5*k22*delta));
       //k4 values
       k41 = Vdot ((vxw+k31*delta), (gamma+k32*delta), (h+k33*delta), CD);
       k42 = Gammadot ((vxw+k31*delta), (gamma+k32*delta), (h+k33*delta), CL);
       k43 = Hdot ((vxw+k31*delta), (gamma+k32*delta));
       k44 = Xdot ((vxw+k31*delta), (gamma+k32*delta));
       //yi+1 for each of the 4 functions to be integrated based on the weighted average of the k values
       y[0] = vxw + (1.0/6.0)*(k11+2.0*k21+2.0*k31+k41)*delta;//integrate Vdot
       y[1] = gamma + (1.0/6.0)*(k12+2.0*k22+2.0*k32+k42)*delta;//integrate Gammadot
       y[2] = h + (1.0/6.0)*(k13+2.0*k23+2.0*k33+k43)*delta;//integrate Hdot
       y[3] = x + (1.0/6.0)*(k14+2.0*k24+2.0*k34+k44)*delta;//integrate Xdot
       y[4] = Ground (y[3]);//plug Xdot into the ground elevation function to calculate the elevation at current X location
       //the y[] array returns the values of the integrations from the void RK function because arrays are treated as pointers
       if (y[3]>6233.58){//ground is level at far end of valley
          y[4]=7480.58;
                         }
       y[5] = Veas (y[0], y[2]);
       }
