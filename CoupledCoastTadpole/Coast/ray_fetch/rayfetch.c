/*
 * =============================================================
 * rayfetch.c - C MEX file to compute fetch polygons
 *
 *
 * This is a MEX-file for MATLAB.
 * Copyright (c) 2020 Taylor Perron
 * =============================================================
 */

#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <stdlib.h>

#define PI 3.1415926535
#define MAXLOOP 100000

double IsBetween(const double astart, const double aend, const double aray) {
    
    // This tests if the ray vector is between the vectors defined by the 
    // two shoreline segments. I think this function wants all three angles
    // to be on the interval [0, 2*pi]
        
    double a1, a2, ar;
    
    // Start by subtracting the start angle from all angles.
    a1 = astart - astart;
    a2 = aend - astart;
    ar = aray - astart;

    // If any angles are now negative, add 360deg to them.
    if (a2 < 0) {
        a2 = a2 + 2*PI;
    } else if (a2 == 0) { // special case: if the prev and next angles are the same, the shoreline doubles back, so make the next angle 360
        a2 = 2*PI;
    }
    
    if (ar < 0) {
        ar = ar + 2*PI;
    }
    
    // Now check to see if the ray angle is smaller than the end angle.
    if (ar <= a2) {
        return 1;
    } else {
        return 0;
    }
    
} // end IsBetween


double IsNotBehind(const int p, const int n, const int k, const double thetaprev, const double thetanext, const double thetaray, double CoastAz[]) {

    // Tests if ray points in front of (or parallel to) the coast, and not behind the coast

    double thetacoast, thetapn, thetanp, dx, dy;
    double irel[] = {0,1,1,1,0,-1,-1,-1};
    double jrel[] = {1,1,0,-1,-1,-1,0,1};

//     // make all the input angles fall in the range [0, 2*pi] (this should already be true for thetaray)
//     if (thetaprev < 0) {thetaprev = thetaprev + 2*PI;} else if (thetaprev > 2*PI) {thetaprev = thetaprev - 2*PI;}
//     if (thetanext < 0) {thetanext = thetanext + 2*PI;} else if (thetanext > 2*PI) {thetanext = thetanext - 2*PI;}

    
    if (p == n) { // special case: if the previous and next shoreline points are the same
        // The azimuth of the normal vector to the coast (pointing into the water)
        // points in the direction opposite from the angle from i,j to the previous (or next) point
        thetacoast = thetaprev + PI;

        // min and max azimuths of rays that point into the water or parallel to the coast
        thetapn = thetacoast - 0.5*PI;
        thetanp = thetacoast + 0.5*PI;
                
    } else { // general case: previous and next shoreline points are different
        // get azimuth from previous shoreline point to next shoreline point
        // because the coast is ordered CW (if i is positive down), the coast normal vector is +90 degrees from this azimuth
//         dx = jrel[n] - jrel[p];
//         dy = irel[n] - irel[p];
        
        // That didn't work...the coast vector pointed towards the land, 180 deg from the correct direction.
        // Do we want next to previous?
        dx = jrel[p] - jrel[n];
        dy = irel[p] - irel[n];
        thetapn = atan2(dy,dx); 
        
        thetanp = thetapn + PI; // azimuth from next shoreline point to previous shoreline point
        thetacoast = thetapn + 0.5*PI;
    }
    
//     // make the azimuth of the normal vector to the coast fall in the range [-pi, pi]
//     if (thetacoast < -PI) {thetacoast = thetacoast + 2*PI;} else if (thetacoast > PI) {thetacoast = thetacoast - 2*PI;}

    // make the azimuth of the normal vector to the coast fall in the range [0, 2*pi]
    if (thetacoast < 0) {thetacoast = thetacoast + 2*PI;} else if (thetacoast > 2*PI) {thetacoast = thetacoast - 2*PI;}
    
    // and save it
    CoastAz[8*p+n] = thetacoast;

    // make all the input angles fall in the range [0, 2*pi] (this should already be true for thetaray)
    if (thetapn < 0) {thetapn = thetapn + 2*PI;} else if (thetapn > 2*PI) {thetapn = thetapn - 2*PI;}
    if (thetanp < 0) {thetanp = thetanp + 2*PI;} else if (thetanp > 2*PI) {thetanp = thetanp - 2*PI;}
    
    // mexPrintf("p=%d n=%d thetacoast=%f k=%d thetapn=%f thetanp=%f thetaray=%f IsBetween=%f\n",p,n,thetacoast,k,thetapn,thetanp,thetaray,IsBetween(thetapn,thetanp,thetaray));
    
    return IsBetween(thetapn,thetanp,thetaray); // If the ray points into the water or parallel to the coast, return 1; otherwise, return 0
    
} // end IsNotBehind


void GetDiDj(double Di[] , double Dj[], double RayAz[], double CoastAz[], const int nrays, const double delta, double isRayAllowed[]) {
    
    int k, p, n;
    double theta, dtheta, thetaprev, thetanext; 
    double angle[] = {0,0.25*PI,0.5*PI,0.75*PI,PI,1.25*PI,1.5*PI,1.75*PI,2*PI};
            
    dtheta = 2*PI/nrays;
    theta = 0.5*dtheta; // 0;
    
    for (k=0; k<nrays; k++) {
        
        Di[k] = delta*sin(theta);
        Dj[k] = delta*cos(theta);

        RayAz[k] = theta;
//         // make the ray azimuths in the return argument fall in the range [-pi, pi]
//         if (theta < -PI) {RayAz[k] = theta + 2*PI;} else if (theta > PI) {RayAz[k] = theta - 2*PI;}

//         mexPrintf("%f %f %f\n",Di[k],Dj[k],theta-dtheta);
        
        // now let's make an array telling us whether each ray direction is permitted
        // given each possible combination of adjoining shoreline segments.
        // Shoreline segments can be (with numbers indicating indices):
        //  5  6  7
        //   \ | /
        // 4 -   - 0
        //   / | \
        //  3  2  1

        for (p=0; p<8; p++) { // prev shoreline segment
            thetaprev = angle[p];
            for (n=0; n<8; n++) { // next shoreline segment
                thetanext = angle[n];
                
                // is ray k between these two shoreline segments? 

                // Indexing into 3D array: M[layer*nrows*ncols + nrows*j + i]
//                 isRayAllowed[k*8*8 + 8*p + n] = IsBetween(thetaprev,thetanext,theta);
                isRayAllowed[k*8*8 + 8*p + n] = IsBetween(thetaprev,thetanext,theta) * IsNotBehind(p,n,k,thetaprev,thetanext,theta,CoastAz);

//                 mexPrintf("p=%d n=%d k=%d thetap=%03.0f thetan=%03.0f theta=%03.0f isBetween=%g isRayAllowed=%g\n",p,n,k,thetaprev*180/PI,thetanext*180/PI,theta*180/PI,isBetween,isRayAllowed[k*8*8 + 8*p + n]);            
            }
        }
    theta = theta + dtheta;   
    }
        
} // end GetDiDj()


double CCW(const double xa, const double ya, const double xb, const double yb, const double xc, const double yc) 
{
    if ( (yc-ya)*(xb-xa) > (yb-ya)*(xc-xa) ) {
        return 1;
    } else {
        return 0;
    }
} // end CCW()

double DoIntersect(const double xa, const double ya, const double xb, const double yb, const double xc, const double yc, const double xd, const double yd) 
{ 
    // def ccw(A,B,C):
    // 	return (C.y-A.y)*(B.x-A.x) > (B.y-A.y)*(C.x-A.x)
    // 
    // def intersect(A,B,C,D):
    // 	return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)
    // NOTE THAT THIS DOES NOT HANDLE COLINEARITY, SO IF ONE ENDPOINT FALLS ON 
    // THE OTHER SEGMENT, OR IF THE TWO SETS OF POINTS ARE THE SAME, WE WILL MISS IT!

    if ( CCW(xa,ya,xc,yc,xd,yd) != CCW(xb,yb,xc,yc,xd,yd) && CCW(xa,ya,xb,yb,xc,yc) != CCW(xa,ya,xb,yb,xd,yd) ) {
        return 1;
    } else {
        return 0;
    }

//     // Find the four orientations needed for general and 
//     // special cases 
//     int o1 = orientation(p1, q1, p2); 
//     int o2 = orientation(p1, q1, q2); 
//     int o3 = orientation(p2, q2, p1); 
//     int o4 = orientation(p2, q2, q1); 
//   
//     // General case 
//     if (o1 != o2 && o3 != o4) 
//         return 1; 
//   
//     // Special Cases 
//     // p1, q1 and p2 are colinear and p2 lies on segment p1q1 
//     if (o1 == 0 && onSegment(p1, p2, q1)) return true; 
//   
//     // p1, q1 and q2 are colinear and q2 lies on segment p1q1 
//     if (o2 == 0 && onSegment(p1, q2, q1)) return true; 
//   
//     // p2, q2 and p1 are colinear and p1 lies on segment p2q2 
//     if (o3 == 0 && onSegment(p2, p1, q2)) return true; 
//   
//      // p2, q2 and q1 are colinear and q1 lies on segment p2q2 
//     if (o4 == 0 && onSegment(p2, q1, q2)) return true; 
//   
//     return 0; // Doesn't fall in any of the above cases 

} // end DoIntersect()



// void GetFetchPolygon(size_t K, size_t J, size_t ns, const int c, double fpi[], double fpj[], double coastnormal[], double Di[], double Dj[], double CoastAz[], double M[], double si[], double sj[], const int nrays, const double delta, double isRayAllowed[], double ShorelineIdx[]) {
void GetFetchPolygon(size_t K, size_t J, size_t ns, const int c, double fpi[], double fpj[], double coastnormal[], double Di[], double Dj[], double CoastAz[], double M[], double si[], double sj[], const int nrays, const double delta, double isRayAllowed[], double Sidx[], double siprev[], double sinext[], double sjprev[], double sjnext[]) {
// This version uses the following algorithm:
//     in each step, calculate the fractional coordinate we'll land at next
//     find the 4 nearest grid points.
//     If none of the points is land, we're in water and won't cross a shoreline; advance
//     If all 4 of those points are land, we must have crossed the shoreline; don't advance and terminate the ray at the point we were already at
//     If 1, 2, or 3 of the points are land, check if we'll intersect a shoreline segment
    
//     int Muu, Mdd, Mud, Mdu, cnext, cprev, cM, cMprev, cMnext, i, j, iM, jM, iMp, jMp, iMn, jMn, iup, jup, idown, jdown, in, ip, jn, jp, dip, din, djp, djn, p, n, k, landho, loops; // i and j are the point we're on and k is the ray we're on
    int Muu, Mdd, Mud, Mdu, cnext, cprev, cM, i, j, iM, jM, iMp, jMp, iMn, jMn, iup, jup, idown, jdown, in, ip, jn, jp, dip, din, djp, djn, p, n, k, landho, loops; // i and j are the point we're on and k is the ray we're on
    int nidx[] = {5,4,3,6,8,2,7,0,1};
    double id, jd, idn, jdn, dy, dx; // decimal row and column coordinates along the ray

    i = (int) (si[c]);
    j = (int) (sj[c]);
    
//     if (c == ns-1) {cnext = 0;} else {cnext = c+1;} // This would be for a shoreline ordered opposite from Rose's code
//     if (c == 0) {cprev = ns-1;} else {cprev = c-1;}

    if (c == ns-1) {cprev = 0;} else {cprev = c+1;} // Because Rose's code orders shoreline CW if i is positive down, CCW if i is positive up 
    if (c == 0) {cnext = (int)ns-1;} else {cnext = c-1;}

    in = (int) (si[cnext]);
    ip = (int) (si[cprev]);

    jn = (int) (sj[cnext]);
    jp = (int) (sj[cprev]);
    
    // The shoreline is ordered CW (if i increases down) and our rays are ordered CW (if i increases down), so
    // The prev angle is the vector from i,j to ip,jp
    // The next angle is the vector from i,j to in,jn
    // Find which neighbor index each of these segments points to
    //  5  6  7
    //   \ | /
    // 4 -   - 0
    //   / | \
    //  3  2  1
    
    dip = ip - i + 1; // These are the zero-based row (for i) and column (for j) indices
    djp = jp - j + 1;
    din = in - i + 1;
    djn = jn - j + 1;
    
    p = nidx[3*djp + dip];
    n = nidx[3*djn + din];

    // save the azimuth of the coast-normal vector in the output argument
    coastnormal[c] = CoastAz[8*p+n];
    
//     mexPrintf("i=%d j=%d p=%d n=%d\n",i,j,p,n);
    
    for (k=0; k<nrays; k++) { // loop through the ray directions

        // set the starting location
        id = si[c];
        jd = sj[c];

        // Use the orientations of the shoreline segments
        // leading to and from this point to limit which rays we try. 
//         mexPrintf("\tk=%d theta=%03.0f isRayAllowed=%g\n",k,180/PI*atan2(Di[k],Dj[k]),isRayAllowed[k*8*8 + 8*p + n]);
                
        if (isRayAllowed[k*8*8 + 8*p + n]) { // if ray k points inside or parallel to the shoreline
                        
            landho = 0;
            loops = 0;

            // march along the ray
            while (!landho && loops<MAXLOOP) { // keep going as long as we're on the grid and haven't hit new land

                idn = id + Di[k]; // the next location along the ray
                jdn = jd + Dj[k];

                if (idn<1 || idn>K || jdn<1 || jdn>J) { // if the next location is off the grid, landho and don't go there
                    landho = 1;
                } else { // the next location is still on the grid
                    
                    if (loops == 0) { // if we're still at the starting point, just advance. This assumes our step length delta is less than the grid spacing!
                        id = idn; // advance
                        jd = jdn;
                    } else {
                        // find the 4 nearest points and their land or water values
                        // (Note that we won't get 4 unique points if idn,jdn falls
                        // exactly on a grid point. This may lead to missed shoreline 
                        // intersections.)
                        iup = (int) ceil(idn); // find the nearest point by rounding to the nearest integer
                        jup = (int) ceil(jdn);
                        idown = (int) floor(idn); // find the nearest point by rounding to the nearest integer
                        jdown = (int) floor(jdn);

                        Muu = (int) M[K*(jup-1)+(iup-1)];
                        Mdd = (int) M[K*(jdown-1)+(idown-1)];
                        Mud = (int) M[K*(jdown-1)+(iup-1)];
                        Mdu = (int) M[K*(jup-1)+(idown-1)];

                        switch ( Muu + Mdd + Mdu + Mud ) {

                            case 0 : // all 4 points are land
                                landho = 1; // terminate
                                break;

                            case 1 : // 3, 2, or 1 neighbors are land
                            case 2 :
                            case 3 :
                                // Test each neighbor to see if it's land and
                                // on the shoreline. If on the shoreline, find
                                // the shoreline segments leading to and from
                                // the neighbor point and test whether this
                                // step would cross each one. If we find a 
                                // shoreline segment that would be crossed, 
                                // terminate and don't bother checking any others.
                                // If we find none that would be crossed,
                                // advance.

                                if (!Muu) { // If Muu is land
                                    iM = iup;
                                    jM = jup;
//                                     cM = (int) ShorelineIdx[K*(jM-1)+(iM-1)]; // this is the 1-based index into si, sj
                                    cM = (int) Sidx[K*(jM-1)+(iM-1)]; // this is the 1-based index into sinext, siprev, sjnext, sjprev
                                    if (cM) { // and if it's on the shoreline
                                        // Find the 2 points that define its 
                                        // previous and next shoreline segments
                                        cM = cM-1; // Switch to zero-based index
//                                         if (cM == ns-1) {cMprev = 0;} else {cMprev = cM+1;} // Because Rose's code orders shoreline CW if i is positive down, CCW if i is positive up 
//                                         if (cM == 0) {cMnext = (int)ns-1;} else {cMnext = cM-1;}

//                                         iMn = (int) (si[cMnext]);
//                                         iMp = (int) (si[cMprev]);
// 
//                                         jMn = (int) (sj[cMnext]);
//                                         jMp = (int) (sj[cMprev]);

                                        iMn = (int) (sinext[cM]);
                                        iMp = (int) (siprev[cM]);

                                        jMn = (int) (sjnext[cM]);
                                        jMp = (int) (sjprev[cM]);
                                        
                                        // If the ray would intersect either of 
                                        // those segments
                                        if ( DoIntersect(id,jd,idn,jdn,(double)iM,(double)jM,(double)iMp,(double)jMp) ) { // if ray intersects prev segment
                                            landho = 1;
                                            break;
                                        } else if ( DoIntersect(id,jd,idn,jdn,(double)iM,(double)jM,(double)iMn,(double)jMn) ) { // if ray intersects next segment
                                            landho = 1;
                                            break;
                                        } // If no intersections, continue

                                    } // If it's not on the shoreline, continue
                                } // end checking Muu

                                if (!Mdd) { // Now check Mdd
                                    iM = idown;
                                    jM = jdown;
//                                     cM = (int) ShorelineIdx[K*(jM-1)+(iM-1)]; // this is the 1-based index into si, sj
//                                     if (cM) { // and if it's on the shoreline
//                                         // Find the 2 points that define its 
//                                         // previous and next shoreline segments
//                                         cM = cM-1; // Switch to zero-based index
//                                         if (cM == ns-1) {cMprev = 0;} else {cMprev = cM+1;} // Because Rose's code orders shoreline CW if i is positive down, CCW if i is positive up 
//                                         if (cM == 0) {cMnext = (int)ns-1;} else {cMnext = cM-1;}
// 
//                                         iMn = (int) (si[cMnext]);
//                                         iMp = (int) (si[cMprev]);
// 
//                                         jMn = (int) (sj[cMnext]);
//                                         jMp = (int) (sj[cMprev]);

                                    cM = (int) Sidx[K*(jM-1)+(iM-1)]; // this is the 1-based index into sinext, siprev, sjnext, sjprev
                                    if (cM) { // and if it's on the shoreline
                                        // Find the 2 points that define its 
                                        // previous and next shoreline segments
                                        cM = cM-1; // Switch to zero-based index

                                        iMn = (int) (sinext[cM]);
                                        iMp = (int) (siprev[cM]);

                                        jMn = (int) (sjnext[cM]);
                                        jMp = (int) (sjprev[cM]);
                                                                        
                                        // If the ray would intersect either of 
                                        // those segments
                                        if ( DoIntersect(id,jd,idn,jdn,(double)iM,(double)jM,(double)iMp,(double)jMp) ) { // if ray intersects prev segment
                                            landho = 1;
                                            break;
                                        } else if ( DoIntersect(id,jd,idn,jdn,(double)iM,(double)jM,(double)iMn,(double)jMn) ) { // if ray intersects next segment
                                            landho = 1;
                                            break;
                                        } // If no intersections, continue

                                    } // If it's not on the shoreline, continue                                
                                }

                                if (!Mud) { // Now check Mud
                                    iM = iup;
                                    jM = jdown;
//                                     cM = (int) ShorelineIdx[K*(jM-1)+(iM-1)]; // this is the 1-based index into si, sj
//                                     if (cM) { // and if it's on the shoreline
//                                         // Find the 2 points that define its 
//                                         // previous and next shoreline segments
//                                         cM = cM-1; // Switch to zero-based index
//                                         if (cM == ns-1) {cMprev = 0;} else {cMprev = cM+1;} // Because Rose's code orders shoreline CW if i is positive down, CCW if i is positive up 
//                                         if (cM == 0) {cMnext = (int)ns-1;} else {cMnext = cM-1;}
// 
//                                         iMn = (int) (si[cMnext]);
//                                         iMp = (int) (si[cMprev]);
// 
//                                         jMn = (int) (sj[cMnext]);
//                                         jMp = (int) (sj[cMprev]);

                                    cM = (int) Sidx[K*(jM-1)+(iM-1)]; // this is the 1-based index into sinext, siprev, sjnext, sjprev
                                    if (cM) { // and if it's on the shoreline
                                        // Find the 2 points that define its 
                                        // previous and next shoreline segments
                                        cM = cM-1; // Switch to zero-based index

                                        iMn = (int) (sinext[cM]);
                                        iMp = (int) (siprev[cM]);

                                        jMn = (int) (sjnext[cM]);
                                        jMp = (int) (sjprev[cM]);
                                    
                                        // If the ray would intersect either of 
                                        // those segments
                                        if ( DoIntersect(id,jd,idn,jdn,(double)iM,(double)jM,(double)iMp,(double)jMp) ) { // if ray intersects prev segment
                                            landho = 1;
                                            break;
                                        } else if ( DoIntersect(id,jd,idn,jdn,(double)iM,(double)jM,(double)iMn,(double)jMn) ) { // if ray intersects next segment
                                            landho = 1;
                                            break;
                                        } // If no intersections, continue

                                    } // If it's not on the shoreline, continue                                
                                }

                                if (!Mdu) { // Now check Mdu
                                    iM = idown;
                                    jM = jup;
//                                     cM = (int) ShorelineIdx[K*(jM-1)+(iM-1)]; // this is the 1-based index into si, sj
//                                     if (cM) { // and if it's on the shoreline
//                                         // Find the 2 points that define its 
//                                         // previous and next shoreline segments
//                                         cM = cM-1; // Switch to zero-based index
//                                         if (cM == ns-1) {cMprev = 0;} else {cMprev = cM+1;} // Because Rose's code orders shoreline CW if i is positive down, CCW if i is positive up 
//                                         if (cM == 0) {cMnext = (int)ns-1;} else {cMnext = cM-1;}
// 
//                                         iMn = (int) (si[cMnext]);
//                                         iMp = (int) (si[cMprev]);
// 
//                                         jMn = (int) (sj[cMnext]);
//                                         jMp = (int) (sj[cMprev]);

                                    cM = (int) Sidx[K*(jM-1)+(iM-1)]; // this is the 1-based index into sinext, siprev, sjnext, sjprev
                                    if (cM) { // and if it's on the shoreline
                                        // Find the 2 points that define its 
                                        // previous and next shoreline segments
                                        cM = cM-1; // Switch to zero-based index

                                        iMn = (int) (sinext[cM]);
                                        iMp = (int) (siprev[cM]);

                                        jMn = (int) (sjnext[cM]);
                                        jMp = (int) (sjprev[cM]);

                                    // If the ray would intersect either of 
                                        // those segments
                                        if ( DoIntersect(id,jd,idn,jdn,(double)iM,(double)jM,(double)iMp,(double)jMp) ) { // if ray intersects prev segment
                                            landho = 1;
                                            break;
                                        } else if ( DoIntersect(id,jd,idn,jdn,(double)iM,(double)jM,(double)iMn,(double)jMn) ) { // if ray intersects next segment
                                            landho = 1;
                                            break;
                                        } // If no intersections, continue

                                    } // If it's not on the shoreline, continue                                
                                }

                                // If we've tried all the neighbors and found
                                // no intersections with shoreline segments, 
                                // advance.
                                id = idn;
                                jd = jdn;
                                break;

                            case 4 : // all 4 points are water
                                id = idn; // advance
                                jd = jdn;
                                break;

                        } // end switch                
                    } // end checking if we're still at the starting point
                } // end checking that the next location was on the grid
                loops++; // count how long we've been in the while loop
            } // end while loop

        } // if ray k points outside the shoreline, leave the starting point as the ray endpoint


        // Save the ray endpoint (which could be the same as the starting point)
        // as the fetch polygon vertex for this direction
        // Note that these fetches might tend to be a bit short because the ray terminates
        // before it crosses onto land, but the angles will be true
        fpi[ns*k + c] = id;
        fpj[ns*k + c] = jd;

    } // done with this ray

} // end GetFetchPolygon()



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    
//     mxArray *Diptr, *Djptr, *isRayAllowedptr, *ShorelineIdxptr, *CoastAzptr;
    mxArray *Diptr, *Djptr, *isRayAllowedptr, *CoastAzptr;
    double *M, *fpi, *fpj, *coastnormal, *si, *sj, *Sidx, *siprev, *sinext, *sjprev, *sjnext, *Di, *Dj, *RayAz, *CoastAz, *isRayAllowed, *ShorelineIdx, *nraysptr, *deltaptr, delta;
    int c, i, j, nrays, ndims=3;
    size_t K, J, ns;
    mwSize dims[]={8, 8, 0};
    
    // Get pointers to inputs
    M = (double *)mxGetPr(prhs[0]); // lake matrix (water = 1, land = 0) 
    si = (double *)mxGetPr(prhs[1]); // row indices of points on shoreline, ordered CCW, first point != last point
    sj = (double *)mxGetPr(prhs[2]); // column indices of points on shoreline, ordered CCW, first point != last point
    nraysptr = (double *)mxGetPr(prhs[3]); /* number of rays around the circle */
    deltaptr = (double *)mxGetPr(prhs[4]); // distance increment along each ray, in cells
    Sidx = (double *)mxGetPr(prhs[5]); // 1-based indices of each shoreline point in Sprev and Snext
    siprev = (double *)mxGetPr(prhs[6]); // vector of 1-based row indices of previous points in the lake matrix
    sinext = (double *)mxGetPr(prhs[7]); // vector of 1-based row indices of next points in the lake matrix
    sjprev = (double *)mxGetPr(prhs[8]); // vector of 1-based row indices of previous points in the lake matrix
    sjnext = (double *)mxGetPr(prhs[9]); // vector of 1-based row indices of next points in the lake matrix
    
    nrays = (int) nraysptr[0];
    delta = deltaptr[0];
    dims[2]= nrays;    
    
    // Get dimensions of input arrays
    K=mxGetM(prhs[0]);
    J=mxGetN(prhs[0]);
    ns=mxGetM(prhs[1]);
    
    // Create arrays for return arguments
    fpi = (double *)mxGetPr(plhs[0]= mxCreateDoubleMatrix(ns, nrays, mxREAL)); // row indices of fetch polygon vertices
    fpj = (double *)mxGetPr(plhs[1]= mxCreateDoubleMatrix(ns, nrays, mxREAL)); // column indices of fetch polygon vertices
    RayAz = (double *)mxGetPr(plhs[2]= mxCreateDoubleMatrix(1, nrays, mxREAL)); // azimuths of rays, in radians. CCW from east if i is positive up, CW from east if i is positive down
    coastnormal = (double *)mxGetPr(plhs[3]= mxCreateDoubleMatrix(ns, 1, mxREAL)); // column vector of azimuths of coast-normal vectors, in radians. CCW from east if i is positive up, CW from east if i is positive down
    
    // Create internally used arrays
    Di = (double *)mxGetPr(Diptr= mxCreateDoubleMatrix(1, nrays, mxREAL)); /* vertical distance increment along a ray */
    Dj = (double *)mxGetPr(Djptr= mxCreateDoubleMatrix(1, nrays, mxREAL)); /* horizontal distance increment along a ray */
    isRayAllowed = (double *)mxGetPr(isRayAllowedptr= mxCreateNumericArray(ndims, dims, mxDOUBLE_CLASS, mxREAL));
    // ShorelineIdx = (double *)mxGetPr(ShorelineIdxptr= mxCreateDoubleMatrix(K, J, mxREAL));
    CoastAz = (double *)mxGetPr(CoastAzptr= mxCreateDoubleMatrix(8, 8, mxREAL)); /* coast azimuths for any combo of neighbor points */
    
    // Find the horizontal and vertical distance increments in each ray orientation
    GetDiDj(Di,Dj,RayAz,CoastAz,nrays,delta,isRayAllowed);

    // Populate ShorelineIdx, which is zero where points are not on the 
    // shoreline and contains the 1-based index into si, sj where points 
    // are on the shoreline
//     for (c=0; c<ns; c++) {            
//         ShorelineIdx[K*((int)sj[c]-1)+((int)si[c]-1)] = c+1;
//     }
    
    // Get fetch polygon vertices for each shoreline point
    for (c=0; c<ns; c++) {
        
//         GetFetchPolygon(K, J, ns, c, fpi, fpj, coastnormal, Di, Dj, CoastAz, M, si, sj, nrays, delta, isRayAllowed,ShorelineIdx);
        GetFetchPolygon(K, J, ns, c, fpi, fpj, coastnormal, Di, Dj, CoastAz, M, si, sj, nrays, delta, isRayAllowed, Sidx, siprev, sinext, sjprev, sjnext);

    }
    
    mxDestroyArray(Diptr);
    mxDestroyArray(Djptr);            
    mxDestroyArray(isRayAllowedptr);            
//     mxDestroyArray(ShorelineIdxptr);            
    mxDestroyArray(CoastAzptr);            

} // end mexFunction