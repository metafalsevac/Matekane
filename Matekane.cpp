/* Final Assignment - Matekane Airstrip
**created December 8, 2011 by metafalsevac
**
*/
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

#include <limits>

// Physical constants. These used to be single-letter #define macros (a, R, W, m, S, g),
// which silently rewrite every later identifier with the same name; typed constants are safe.
constexpr double T0 = 518.67;     //sea level temperature in Rankine
constexpr double P0 = 2116.224;   //sea level pressure in lb/ft^2
constexpr double a = -0.003566;   //Temperature Lapse Rate in Rankine/ft
constexpr double R = 1716.0;      //Specific Gas Constant for air in ft*lbf/slug/degree Rankine
constexpr double W = 2712.5;      //modified weight in pounds (3100 lb reduced by 12.5%)
constexpr double g = 32.2;        //gravitational acceleration in ft/sec^2
constexpr double m = W / g;       //modified mass in slugs (= 84.2391...)
constexpr double S = 195.75;      //modified reference area in sq ft (174 sq ft +12.5%)
constexpr double CLMax = 2.195;   //maximum (stall) lift coefficient
constexpr double T0SL = 714.9375; //modified max sea level thrust in pounds (635.5 lb +12.5%)
constexpr double aSL = 0.0334;    //airspeed variation coefficient for max thrust at sea level in lbs*sec^2/ft^2
constexpr double Vxw0 = 83.0;     //initial condition for Vxw (true airspeed) in ft/s
constexpr double X0 = 0.0;        //initial x distance
constexpr double gamma0 = 0.0;    //initial flight path angle gamma (radians)
constexpr double H0 = 7649.63;    //initial altitude in feet (= Ground(0))
constexpr double XLevel = 6233.58;  //beyond this x the ground is level
constexpr double HLevel = 7480.58;  //elevation of the level ground
constexpr double XEnd = 10000.0;    //integrate until this x
constexpr double TMax = 600.0;      //safety limit on simulated time (s) to avoid infinite loops
constexpr double DT = 0.1;          //RK4 time step (s)

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
                                                                                           //this function is void because it needs to return 6 values
                                                                                           //it does so using arrays, which are treated as pointers
int main(){
double y[6], CLmin, CLmax, step;//array used to return results of integration
double CL, CD;

ofstream outfile ("output.txt");//open output file
outfile << "# x(ft)  h(ft)  V(ft/s)  ground(ft)  gamma(rad)  Veas(ft/s)  CL" << endl;

cout << "Takeoff From Matekane Airstrip\n\n";
cout << "This program will simulate takeoff from the Matekane Airstrip \nat various lift coefficients.\n\nInput minimum lift coefficient: ";
cin >> CLmin;
cout << "\nInput maximum lift coefficient: ";
cin >> CLmax;
cout << "\nInput lift coefficient step size: ";
cin >> step;
cout << endl;

//input validation: a zero/negative step used to loop forever
if (!cin || !(step > 0.0) || CLmax < CLmin) {
   cout << "Invalid input: need CLmin <= CLmax and step > 0." << endl;
   return 1;
   }

//run RK for multiple CL values. An integer counter is used because accumulating
//"CL += step" in floating point skips the last value (e.g. 1.0..2.0 step 0.1 never ran CL = 2.0)
int nCL = (int)floor((CLmax - CLmin)/step + 1e-9) + 1;
for (int i = 0; i < nCL; ++i) {
CL = CLmin + i*step;
if (CL > CLMax + 1e-12) {
   cout << "CL = " << CL << " exceeds CLMax = " << CLMax << " (stalled wing) - skipped" << endl << endl;
   continue;
   }
y[0] = Vxw0;//must reset ALL initial conditions each iteration
y[1] = gamma0;
y[2] = H0;
y[3] = X0;
y[4] = Ground(X0);
y[5] = Veas(Vxw0, H0);//previously not reset, so each run's first output row carried the previous run's final Veas
CD = 0.025 + 0.056*CL*CL;
double t = 0.0;
bool crashed = false, diverged = false;

    while (y[3]<=XEnd) {//while still in the valley, the while loop actually runs the integration past the far end of the valley to determine if the plane continues to gain altitude
      outfile << y[3] << "     " << y[2] << "     " << y[0] << "      " << y[4] << "     " << y[1] << "     " << y[5] << "     " << CL << endl;
      double xPrev = y[3], hPrev = y[2], gPrev = y[4];

      RK (y[0], y[1], y[2], y[3], y, CD, CL);//call runge kutta function
      t += DT;

      if (!isfinite(y[0]) || !isfinite(y[1]) || y[0] <= 0.0 || t > TMax) {//airspeed lost / time limit: avoid NaN or an endless loop
         cout << "CL = " << CL << " NO RESULT (airspeed lost or time limit reached) at t = " << t << " s" << endl << endl;
         diverged = true;
         break;
         }

      if ((y[4]>=y[2])){//did plane crash?
         //the crash is only detected after the step that went below ground; interpolate the touchdown point
         double c0 = hPrev - gPrev, c1 = y[2] - y[4];
         double xTouch = xPrev + (y[3] - xPrev) * (c0 / (c0 - c1));
         cout << "CL = " << CL << " CRASH!!!" << endl;
         cout << "Vxw = " << y[0] << " gamma = " << y[1] << " Alt = " << y[2] << " x = " << y[3] << " Elev = " << y[4] << "\nEquiv Airspeed = " << y[5] << "\nTouchdown (interpolated) x = " << xTouch << endl << endl;//coordinates and flight conditions at crash
         crashed = true;
         break;
                      }//if

                           }//while
      if (!crashed && !diverged && (y[3]>XLevel)&&y[2]>y[4]){
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
       if (x > XLevel) return HLevel;//ground is level at far end of valley (polynomial is continuous there: 7480.579)
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
       
double Hdot (double vxw, double gamma) {//vertical speed
       double hdot = vxw*sin(gamma);
       return hdot;
       }
       
double Xdot (double vxw, double gamma) {//horizontal speed
       double xdot = vxw*cos(gamma);
       return xdot;
       }
       
void RK (double vxw, double gamma, double h, double x, double y[], double CD, double CL) {//runge kutta. this function is void and uses Arrays to return results of integration
       double delta = DT;//increment
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
       y[5] = Veas (y[0], y[2]);
       }
