// Axisymmetric (r,z) Boussinesq convection prototype.
// A buoyancy "bump" imposed at the bottom boundary drives a vorticity/velocity
// field through a relaxation (damped Gauss-Seidel) solver. The vorticity boundary
// condition is not solved numerically but taken from an analytic plume-scaling
// law (Om ~ nu^-1/3 * B^2/3) -- hence "AnalyticConvection". See README.md.

#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <math.h>
#include <time.h>

// ---- grid & physical constants -------------------------------------------

constexpr int GRID_NR = 100; // number of radial grid points
constexpr int GRID_NZ = 100; // number of vertical grid points

constexpr double CELL_SIZE_R = 0.1; // radial cell size, meters
constexpr double CELL_SIZE_Z = 0.1; // vertical cell size, meters

// bottom-boundary buoyancy bump: B(r,0) = BUOYANCY_BOUNDARY_PEAK / (|i - GRID_NR/2| + 1)
constexpr double BUOYANCY_BOUNDARY_PEAK = 100.0;
// analytic vorticity-generation scaling law: Om(r,0) = VORTICITY_SCALING_COEFF * nu^(-1/3) * B(r,0)^(2/3)
constexpr double VORTICITY_SCALING_COEFF = 0.02;

constexpr double RELAXATION_RATE = 0.00005; // under-relaxation factor for the pseudo-time solver
constexpr int TIME_STEPS_PER_FRAME = 50;    // relaxation sub-steps executed per rendered frame

constexpr int VECTOR_GRID_SAMPLES = 25; // approx. number of velocity-vector arrows per axis
constexpr double DOMAIN_WIDTH_R = GRID_NR * CELL_SIZE_R; // radial extent of the domain, meters
constexpr double MAX_VECTOR_LENGTH = DOMAIN_WIDTH_R / 30.0; // length of the fastest vector, meters
constexpr int WINDOW_SIZE = 800;        // window width/height, pixels

// ---- simulation state ------------------------------------------------------

double buoyancy[GRID_NR][GRID_NZ];  // B: buoyancy perturbation, m/s^2
double vorticity[GRID_NR][GRID_NZ]; // Om: azimuthal vorticity, 1/s
double velocityZ[GRID_NR][GRID_NZ]; // Uz: vertical velocity, m/s
double velocityR[GRID_NR][GRID_NZ]; // Ur: radial velocity, m/s

double vorticityRHS[GRID_NR][GRID_NZ];
double velocityZRHS[GRID_NR][GRID_NZ];
double buoyancyRHS[GRID_NR][GRID_NZ];

double viscosity;           // nu
double buoyancyDiffusivity; // k
double buoyancyFrequency;   // N (Brunt-Vaisala-like frequency)

// ---- view / UI state --------------------------------------------------------

enum FieldType
{
    FIELD_BUOYANCY,
    FIELD_VORTICITY,
    FIELD_VELOCITY_Z,
    FIELD_VELOCITY_R
};

FieldType activeField = FIELD_BUOYANCY;

double zoomScale = 0.5;
double panX = 5.0;
double panY = 5.0;
double colorScale = 1.0;

bool isRunning = false;
bool showVectors = false;
double vectorLengthScale = 1.0; // user-adjustable multiplier, keys '[' and ']'

// ---- solver -----------------------------------------------------------------

// Relaxes `field` along z for one radial row `i`, iterating the discretized
// 1D equation  -field[j-1] + 2*field[j] - field[j+1] = rhs[j]  with under-relaxation
// factor `relaxFactor`. field[i][0] (bottom boundary) is never touched here, so it
// keeps whatever fixed (Dirichlet) value initSimulation() gave it. The top boundary
// field[i][GRID_NZ-1] gets a zero-gradient (Neumann) condition.
void relaxImplicitZ(int i, double field[GRID_NR][GRID_NZ], int iterations,
                     const double rhs[GRID_NR][GRID_NZ], double relaxFactor)
{
    const double diagonalCoeff = 2.0;
    const double offDiagonalCoeff = 1.0;

    for (int n = 0; n < iterations; n++)
    {
        for (int j = 1; j < GRID_NZ - 1; j++)
        {
            field[i][j] = field[i][j] * (1.0 - relaxFactor)
                        + relaxFactor * (offDiagonalCoeff * (field[i][j + 1] + field[i][j - 1]) + rhs[i][j]) / diagonalCoeff;
        }
        field[i][GRID_NZ - 1] = field[i][GRID_NZ - 2];
    }
}

// Recovers the radial velocity Ur from Uz by integrating the incompressibility
// (continuity) equation in cylindrical coordinates: (1/r) d(r*Ur)/dr + dUz/dz = 0
void computeRadialVelocityFromContinuity()
{
    for (int j = 1; j < GRID_NZ - 1; j++)
    {
        velocityR[0][j] = 0.0; // r = 0 is the symmetry axis
        for (int i = 1; i < GRID_NR - 1; i++)
        {
            velocityR[i][j] = velocityR[i - 1][j]
                - 0.5 * CELL_SIZE_R * (
                    (i * CELL_SIZE_R) * (velocityZ[i][j + 1] - velocityZ[i][j - 1]) / (2.0 * CELL_SIZE_Z) +
                    ((i - 1) * CELL_SIZE_R) * (velocityZ[i - 1][j + 1] - velocityZ[i - 1][j - 1]) / (2.0 * CELL_SIZE_Z));
        }
        for (int i = 1; i < GRID_NR - 1; i++)
            velocityR[i][j] /= (i * CELL_SIZE_R);
    }
}

// Advances the model by one relaxation step: sweeps radial rows from the
// symmetry axis outward, solving (per row, implicitly in z) vorticity transport,
// then the vorticity -> velocity relation for Uz, then buoyancy transport.
// Radial coupling between rows is explicit / Gauss-Seidel (uses the row i-1
// already updated earlier in the same sweep).
void advanceTimestep(int iterations, double relaxFactor)
{
    for (int i = 0; i < GRID_NR; i++)
    {
        if (i == 0)
        {
            // symmetry axis: no radial advection or baroclinic (buoyancy-gradient) term
            for (int j = 1; j < GRID_NZ; j++)
                vorticityRHS[i][j] = 0.0;
            relaxImplicitZ(i, vorticity, iterations, vorticityRHS, relaxFactor);

            for (int j = 1; j < GRID_NZ; j++)
                velocityZRHS[i][j] = -vorticity[i][j] / CELL_SIZE_R * CELL_SIZE_Z * CELL_SIZE_Z
                                    - CELL_SIZE_Z * CELL_SIZE_Z * vorticity[i][j] / ((i + 1) * CELL_SIZE_R);
            relaxImplicitZ(i, velocityZ, iterations, velocityZRHS, relaxFactor * 0.1);
        }
        else
        {
            // vorticity transport: advection (Ur*Om/r) + baroclinic torque from radial buoyancy gradient
            for (int j = 1; j < GRID_NZ; j++)
                vorticityRHS[i][j] = -(velocityR[i][j] * vorticity[i][j] / (CELL_SIZE_R * (i + 1)) * CELL_SIZE_Z * CELL_SIZE_Z / viscosity
                                     + (buoyancy[i][j] - buoyancy[i - 1][j]) / CELL_SIZE_R * CELL_SIZE_Z * CELL_SIZE_Z / viscosity);
            relaxImplicitZ(i, vorticity, iterations, vorticityRHS, relaxFactor);

            // vorticity -> velocity Poisson relation (cylindrical Laplacian of Uz)
            for (int j = 1; j < GRID_NZ; j++)
                velocityZRHS[i][j] = -CELL_SIZE_Z * CELL_SIZE_Z * (
                    (-2.0 * velocityZ[i][j] + velocityZ[i - 1][j] + velocityZ[i + 1][j]) / (CELL_SIZE_R * CELL_SIZE_R)
                    + (vorticity[i][j] - vorticity[i - 1][j]) / CELL_SIZE_R
                    + (vorticity[i][j] + (velocityZ[i][j] - velocityZ[i - 1][j]) / CELL_SIZE_R) / ((i + 1) * CELL_SIZE_R));
            relaxImplicitZ(i, velocityZ, iterations, velocityZRHS, relaxFactor * 0.1);
        }

        // linearized buoyancy transport: diffuses and is forced by vertical advection of the background stratification
        for (int j = 1; j < GRID_NZ; j++)
            buoyancyRHS[i][j] = velocityZ[i][j] * CELL_SIZE_Z * CELL_SIZE_Z * buoyancyFrequency * buoyancyFrequency / buoyancyDiffusivity;
        relaxImplicitZ(i, buoyancy, iterations, buoyancyRHS, relaxFactor);
    }
    computeRadialVelocityFromContinuity();
}

// ---- field metadata ----------------------------------------------------------

double getFieldValue(FieldType field, int i, int j)
{
    switch (field)
    {
    case FIELD_VORTICITY:  return vorticity[i][j];
    case FIELD_VELOCITY_Z: return velocityZ[i][j];
    case FIELD_VELOCITY_R: return velocityR[i][j];
    default:               return buoyancy[i][j];
    }
}

struct FieldInfo
{
    const char* name;
    const char* units;
    const char* colorLegend;
};

FieldInfo getFieldInfo(FieldType field)
{
    switch (field)
    {
    case FIELD_VORTICITY:  return { "Om (vorticity)", "1/s", "color: green = positive, magenta = negative" };
    case FIELD_VELOCITY_Z: return { "Uz (vertical velocity)", "m/s", "color: red = positive, cyan = negative" };
    case FIELD_VELOCITY_R: return { "Ur (radial velocity)", "m/s", "color: red = positive, cyan = negative" };
    default:               return { "B (buoyancy)", "m/s^2", "color: red = positive, cyan = negative" };
    }
}

void getFieldRange(FieldType field, double& minValue, double& maxValue)
{
    minValue = 1e300;
    maxValue = -1e300;
    for (int i = 0; i < GRID_NR; i++)
        for (int j = 0; j < GRID_NZ; j++)
        {
            double value = getFieldValue(field, i, j);
            if (value < minValue) minValue = value;
            if (value > maxValue) maxValue = value;
        }
}

// ---- rendering ----------------------------------------------------------------

void drawBitmapText(float x, float y, void* font, const char* str)
{
    glRasterPos2f(x, y);
    for (const char* c = str; *c != '\0'; c++)
        glutBitmapCharacter(font, *c);
}

// colors positive values warm (red, or green for vorticity) and negative values cool
void setFieldColor(FieldType field, double value)
{
    double v = colorScale * value;
    if (field == FIELD_VORTICITY)
        glColor3f((float)-v, (float)v, (float)-v);
    else
        glColor3f((float)v, (float)-v, (float)-v);
}

void drawField()
{
    for (int i = 0; i < GRID_NR - 1; i++)
    {
        glBegin(GL_TRIANGLE_STRIP);
        for (int j = 0; j < GRID_NZ; j++)
        {
            setFieldColor(activeField, getFieldValue(activeField, i, j));
            glVertex3f((float)(i * CELL_SIZE_R), (float)(j * CELL_SIZE_Z), 0.0f);

            setFieldColor(activeField, getFieldValue(activeField, i + 1, j));
            glVertex3f((float)((i + 1) * CELL_SIZE_R), (float)(j * CELL_SIZE_Z), 0.0f);
        }
        glEnd();
    }
}

// draws r/z axes with 1-meter tick marks and labels (world coordinates already
// are meters: r = i*CELL_SIZE_R, z = j*CELL_SIZE_Z)
void drawAxes()
{
    double maxR = (GRID_NR - 1) * CELL_SIZE_R;
    double maxZ = (GRID_NZ - 1) * CELL_SIZE_Z;
    const double tickStep = 1.0; // meters
    const double tickLen = 0.15;

    glColor3f(1.0, 1.0, 1.0);
    glLineWidth(2.0);
    glBegin(GL_LINES);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f((float)maxR, 0.0, 0.0);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, (float)maxZ, 0.0);
    glEnd();

    glLineWidth(1.0);
    glBegin(GL_LINES);
    for (double rt = 0.0; rt <= maxR + 1e-9; rt += tickStep)
    {
        glVertex3f((float)rt, 0.0f, 0.0f);
        glVertex3f((float)rt, (float)-tickLen, 0.0f);
    }
    for (double zt = 0.0; zt <= maxZ + 1e-9; zt += tickStep)
    {
        glVertex3f(0.0f, (float)zt, 0.0f);
        glVertex3f((float)-tickLen, (float)zt, 0.0f);
    }
    glEnd();

    char buf[32];
    for (double rt = 0.0; rt <= maxR + 1e-9; rt += tickStep)
    {
        sprintf(buf, "%.0f", rt);
        drawBitmapText((float)(rt - 0.08), -0.55f, GLUT_BITMAP_HELVETICA_10, buf);
    }
    for (double zt = 0.0; zt <= maxZ + 1e-9; zt += tickStep)
    {
        sprintf(buf, "%.0f", zt);
        drawBitmapText(-0.85f, (float)(zt - 0.08), GLUT_BITMAP_HELVETICA_10, buf);
    }

    drawBitmapText((float)(maxR + 0.2), -0.3f, GLUT_BITMAP_HELVETICA_12, "r, m");
    drawBitmapText(-0.85f, (float)(maxZ + 0.3), GLUT_BITMAP_HELVETICA_12, "z, m");
}

// draws velocity vectors (Ur,Uz) as small white points with a thin line per
// sample, on a ~VECTOR_GRID_SAMPLES x VECTOR_GRID_SAMPLES grid. Vector length is
// proportional to velocity amplitude: auto-scaled so the fastest sampled vector
// is MAX_VECTOR_LENGTH long (1/30 of the domain width -- raw velocities are far
// smaller than the domain size and would otherwise be invisible), then further
// scaled by the user-adjustable `vectorLengthScale` ('[' / ']'). This flow has a
// highly peaked velocity field (a narrow rising jet above the buoyancy bump vs. a
// near-still surrounding), so a single dominant sample can dwarf the others once
// `vectorLengthScale` is increased; each vector's on-screen length is therefore
// hard-clamped to 90% of the sampling-grid cell size so no arrow can overrun the
// plot regardless of the underlying velocity magnitude or scale setting.
void drawVelocityVectors()
{
    int stepI = GRID_NR / VECTOR_GRID_SAMPLES;
    int stepJ = GRID_NZ / VECTOR_GRID_SAMPLES;
    if (stepI < 1) stepI = 1;
    if (stepJ < 1) stepJ = 1;

    double maxSpeed = 1e-30;
    for (int i = 0; i < GRID_NR; i += stepI)
        for (int j = 0; j < GRID_NZ; j += stepJ)
        {
            double speed = sqrt(velocityR[i][j] * velocityR[i][j] + velocityZ[i][j] * velocityZ[i][j]);
            if (speed > maxSpeed) maxSpeed = speed;
        }

    double vecScale = vectorLengthScale * MAX_VECTOR_LENGTH / maxSpeed;

    double cellSize = (stepI * CELL_SIZE_R < stepJ * CELL_SIZE_Z) ? stepI * CELL_SIZE_R : stepJ * CELL_SIZE_Z;
    double maxArrowLength = 0.9 * cellSize;

    glColor3f(1.0, 1.0, 1.0);
    glPointSize(1.5);
    glBegin(GL_POINTS);
    for (int i = 0; i < GRID_NR; i += stepI)
        for (int j = 0; j < GRID_NZ; j += stepJ)
            glVertex3f((float)(i * CELL_SIZE_R), (float)(j * CELL_SIZE_Z), 0.0f);
    glEnd();

    glLineWidth(1.0);
    glBegin(GL_LINES);
    for (int i = 0; i < GRID_NR; i += stepI)
        for (int j = 0; j < GRID_NZ; j += stepJ)
        {
            double speed = sqrt(velocityR[i][j] * velocityR[i][j] + velocityZ[i][j] * velocityZ[i][j]);
            if (!(speed > 0.0)) continue; // skips exactly-zero and NaN/garbage speeds

            // fmin() picks the finite operand when either side is NaN, so this also
            // guards against occasional numerical spikes in the raw velocity field
            double arrowLength = fmin(speed * vecScale, maxArrowLength);
            double r0 = i * CELL_SIZE_R, z0 = j * CELL_SIZE_Z;
            double r1 = r0 + (velocityR[i][j] / speed) * arrowLength;
            double z1 = z0 + (velocityZ[i][j] / speed) * arrowLength;
            glVertex3f((float)r0, (float)z0, 0.0f);
            glVertex3f((float)r1, (float)z1, 0.0f);
        }
    glEnd();
}

// screen-space overlay: active field name, units, current value range and color key
void drawLegend()
{
    FieldInfo info = getFieldInfo(activeField);
    double minValue, maxValue;
    getFieldRange(activeField, minValue, maxValue);

    int windowWidth = glutGet(GLUT_WINDOW_WIDTH);
    int windowHeight = glutGet(GLUT_WINDOW_HEIGHT);

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0, windowWidth, 0, windowHeight);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    glColor3f(1.0, 1.0, 1.0);
    char buf[160];
    sprintf(buf, "Field: %s", info.name);
    drawBitmapText(10.0f, (float)(windowHeight - 20), GLUT_BITMAP_HELVETICA_12, buf);
    sprintf(buf, "Units: %s", info.units);
    drawBitmapText(10.0f, (float)(windowHeight - 38), GLUT_BITMAP_HELVETICA_12, buf);
    sprintf(buf, "Scale: min=%.4g  max=%.4g", minValue, maxValue);
    drawBitmapText(10.0f, (float)(windowHeight - 56), GLUT_BITMAP_HELVETICA_12, buf);
    drawBitmapText(10.0f, (float)(windowHeight - 74), GLUT_BITMAP_HELVETICA_12, info.colorLegend);
    if (showVectors)
    {
        sprintf(buf, "Vectors: (Ur,Uz), length ~ speed, scale=%.2f ('[' / ']')", vectorLengthScale);
        drawBitmapText(10.0f, (float)(windowHeight - 92), GLUT_BITMAP_HELVETICA_12, buf);
    }

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}

// ---- GLUT callbacks -----------------------------------------------------------

void onDisplay(void)
{
    glClear(GL_COLOR_BUFFER_BIT);
    glLoadIdentity();
    glScalef(zoomScale, zoomScale, zoomScale);
    glTranslatef(panX, panY, 0.0);

    drawField();
    drawAxes();
    if (showVectors)
        drawVelocityVectors();
    drawLegend();

    if (isRunning)
    {
        for (int step = 0; step < TIME_STEPS_PER_FRAME; step++)
            advanceTimestep(1, RELAXATION_RATE);
        glutPostRedisplay();
    }

    glutSwapBuffers();
}

void onKeyboard(unsigned char key, int x, int y)
{
    switch (key)
    {
    case '1': activeField = FIELD_BUOYANCY; break;
    case '2': activeField = FIELD_VORTICITY; break;
    case '3': activeField = FIELD_VELOCITY_Z; break;
    case '4': activeField = FIELD_VELOCITY_R; break;

    case 'q': zoomScale *= 1.01; break;
    case 'e': zoomScale /= 1.01; break;
    case 'w': panY -= 5.0; break;
    case 's': panY += 5.0; break;
    case 'a': panX += 5.0; break;
    case 'd': panX -= 5.0; break;

    case ' ':
        isRunning = !isRunning;
        printf("simulation %s\n", isRunning ? "running" : "paused");
        break;

    case '.': colorScale *= 1.3; printf("colorScale=%e\n", colorScale); break;
    case ',': colorScale /= 1.3; printf("colorScale=%e\n", colorScale); break;

    case 'v': showVectors = !showVectors; break;
    case '[':
        vectorLengthScale /= 1.25;
        if (vectorLengthScale < 0.05) vectorLengthScale = 0.05;
        printf("vectorLengthScale=%.3f\n", vectorLengthScale);
        break;
    case ']':
        vectorLengthScale *= 1.25;
        printf("vectorLengthScale=%.3f\n", vectorLengthScale);
        break;
    }
    glutPostRedisplay();
}

// ---- setup ----------------------------------------------------------------

void initSimulation()
{
    viscosity = 1e-2;
    buoyancyDiffusivity = 1e-2;
    buoyancyFrequency = 0.01;

    for (int i = 0; i < GRID_NR; i++)
    {
        for (int j = 0; j < GRID_NZ; j++)
        {
            buoyancy[i][j] = 0.0;
            vorticity[i][j] = 0.0;
            velocityZ[i][j] = 0.0;
        }

        // bottom boundary (z=0): buoyancy bump peaked at the domain center, decaying
        // radially; the matching vorticity boundary value comes from the analytic
        // plume-scaling law. Row j=0 is never touched by relaxImplicitZ afterwards,
        // so this assignment fixes it as a Dirichlet boundary condition for the run.
        double boundaryBuoyancy = BUOYANCY_BOUNDARY_PEAK / (fabs(i - GRID_NR / 2) + 1);
        buoyancy[i][0] = boundaryBuoyancy;
        vorticity[i][0] = VORTICITY_SCALING_COEFF * pow(viscosity, -1.0 / 3.0) * pow(boundaryBuoyancy, 2.0 / 3.0);
    }

    glClearColor(0.0, 0.1, 0.0, 0.0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1.1, CELL_SIZE_R * GRID_NR, -1.1, CELL_SIZE_Z * GRID_NZ, -10.0, 10.0);
    glMatrixMode(GL_MODELVIEW);
}

int main(int argc, char** argv)
{
    srand((unsigned int)time(NULL));
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(WINDOW_SIZE, WINDOW_SIZE);
    glutInitWindowPosition(0, 0);
    glutCreateWindow("Analytic Convection Solver");
    glutDisplayFunc(onDisplay);
    glutKeyboardFunc(onKeyboard);
    initSimulation();
    glutMainLoop();
    return 0;
}
