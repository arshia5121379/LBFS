#include "udf.h"

/* Constants */
#define E0 226270   /* Maximum electric field in V/m (convert kV/cm to V/m) */
#define K1 6543     /* Electric field constant in x-direction in V/m^2 */
#define K2 13083    /* Electric field constant in y-direction in V/m^2 */
#define CHARGE_DENSITY 1.0e11 /* Charge density in electrons/cm^3 */
#define ELECTRON_CHARGE 1.602e-19 /* Elementary charge in Coulombs */
#define DISCHARGE_TIME 6.7e-5 /* Discharge time in seconds */
#define FREQUENCY 5000        /* Frequency in Hz */
#define X_PLASMA_BOUND 0.001  /* Plasma boundary in x-direction (in meters) */
#define Y_PLASMA_BOUND 0.0015 /* Plasma boundary in y-direction (in meters) */
#define BREAKDOWN_FIELD 3.0e6 /* Breakdown electric field (in V/m) */

/* UDF for source term in momentum equation */
DEFINE_SOURCE(plasma_body_force_x, c, t, dS, eqn)
{
    real x[ND_ND]; /* Array to store cell centroid coordinates */
    real Ex, Ey, E, Fx, Fy;
    real source_x; /* Source term in x-direction */
    
    C_CENTROID(x, c); /* Get the centroid coordinates of the cell */

    /* Calculate the electric field magnitude */
    E = E0 - K1 * x[0] - K2 * x[1];

    /* Check if the electric field is greater than the breakdown threshold */
    if (E > BREAKDOWN_FIELD)
    {
        /* Calculate the electric field components */
        Ex = E * K2 / sqrt(K1 * K1 + K2 * K2);
        Ey = E * K1 / sqrt(K1 * K1 + K2 * K2);

        /* Calculate the body force components */
        Fx = Ex * CHARGE_DENSITY * ELECTRON_CHARGE;
        Fy = Ey * CHARGE_DENSITY * ELECTRON_CHARGE;

        /* Apply the source term only in the x-direction for this function */
        source_x = Fx * DISCHARGE_TIME * FREQUENCY;
    }
    else
    {
        /* If electric field is below the breakdown threshold, no plasma is formed */
        source_x = 0.0;
    }

    dS[eqn] = 0.0; /* Derivative of the source term w.r.t. the variable */
    return source_x; /* Return the source term for the x-direction */
}

DEFINE_SOURCE(plasma_body_force_y, c, t, dS, eqn)
{
    real x[ND_ND]; /* Array to store cell centroid coordinates */
    real Ex, Ey, E, Fx, Fy;
    real source_y; /* Source term in y-direction */

    C_CENTROID(x, c); /* Get the centroid coordinates of the cell */

    /* Calculate the electric field magnitude */
    E = E0 - K1 * x[0] - K2 * x[1];

    /* Check if the electric field is greater than the breakdown threshold */
    if (E > BREAKDOWN_FIELD)
    {
        /* Calculate the electric field components */
        Ex = E * K2 / sqrt(K1 * K1 + K2 * K2);
        Ey = E * K1 / sqrt(K1 * K1 + K2 * K2);

        /* Calculate the body force components */
        Fx = Ex * CHARGE_DENSITY * ELECTRON_CHARGE;
        Fy = Ey * CHARGE_DENSITY * ELECTRON_CHARGE;

        /* Apply the source term only in the y-direction for this function */
        source_y = Fy * DISCHARGE_TIME * FREQUENCY;
    }
    else
    {
        /* If electric field is below the breakdown threshold, no plasma is formed */
        source_y = 0.0;
    }

    dS[eqn] = 0.0; /* Derivative of the source term w.r.t. the variable */
    return source_y; /* Return the source term for the y-direction */
}
