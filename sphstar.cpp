#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <vector>
#include "canvas.h"

typedef struct
{
    double x, y;
    double vx, vy;
    double vx0, vy0;
    double ax, ay;
    double density;
}Particle;

double kernel(double x, double y, double h)
{
    double q = sqrt(x*x + y*y) / h;
    double C = 5 / (14 * M_PI * h*h);
    if (0 <= q && 1 > q) {
        return C * (4 - 6 * q*q + 3 * q*q*q);
    } else if (1 <= q && 2 > q) {
        return C * (8 - 12 * q + 6 * q*q - q*q*q);
    } else {
        return 0;
    }
}

double gradkernel(double x, double y, double h)
{
    double q = sqrt(x*x + y*y) / h;
    double C = 5 / (14 * M_PI * h*h);
    if (0 <= q && 1 > q) {
        return C * (9 * q - 12) * x / (h*h);
    } else if (1 <= q && 2 > q) {
        return C * (12 - 3 * q - 12 / q) * x / (h*h);
    } else {
        return 0;
    }
}

class System
{
    unsigned int Npart;
    unsigned int i, j;
    double dt;
    double h;
    double m;
    double damp;
    double k;
    double R0;
    int npoly;
    double lambda;
    vector<Particle> particles;
    protected:
        void SetParams();
        void Init();
    public:
        System();
        System(unsigned int);
        System(double);
        System(unsigned int, double);
        ~System();
        void Run();
        void Plot(Canvas*);
        void Write(FILE*);
};

System::System()
{
    Npart = 100;
    SetParams();
    Init();
}

System::System(unsigned int np)
{
    Npart = np;
    R0 = Npart / 500.0;
    SetParams();
    Init();
}

System::System(double radius0)
{
    R0 = radius0;
    Npart = 500 * R0;
    SetParams();
    Init();
}

System::System(unsigned int np, double radius0)
{
    Npart = np;
    R0 = radius0;
    SetParams();
    Init();
}

System::~System()
{
}

void System::SetParams()
{
    dt = 0.005;
    h = 0.14;
    m = 1.33333 / Npart;
    damp = 1.0;
    k = 0.25;
    npoly = 1;
    lambda = 1.0;
}

void System::Init()
{
    double radius = 0;
    double phi = 0;

    particles.resize(Npart);

    srand(time(NULL));
    for (i = 0; i < particles.size(); i++) {
        radius = sqrt((double) rand() / RAND_MAX ) * R0;
        phi = (double) rand() / RAND_MAX * 2 * M_PI;
        particles[i].x = cos(phi) * radius;
        particles[i].y = sin(phi) * radius;
        particles[i].vx = 0;
        particles[i].vy = 0;
        particles[i].vx0 = 0;
        particles[i].vy0 = 0;
        particles[i].density = 0;
        particles[i].ax = 0;
        particles[i].ay = 0;
    }
}

void System::Run()
{
    double rho = 0;
    double pressure = 0;
    double delkern = 0;
    double velocity = 0;

    for (i = 0; i < particles.size(); i++) {
        particles[i].density = 0;
        particles[i].ax = 0;
        particles[i].ay = 0;
    }

    for (i = 0; i < particles.size(); i++) {
        particles[i].density += m * kernel(0, 0, h);
        for (j = i + 1; j < particles.size(); j++) {
            rho = m * kernel(
                particles[i].x - particles[j].x, 
                particles[i].y - particles[j].y, 
                h
            );
            particles[i].density += rho;
            particles[j].density += rho;
        }
    }

    for (i = 0; i < particles.size(); i++) {
        particles[i].ax -= damp * particles[i].vx + lambda * particles[i].x;
        particles[i].ay -= damp * particles[i].vy + lambda * particles[i].y;
        for (j = i + 1; j < particles.size(); j++) {
            pressure = -m * (
                k * pow(particles[i].density, (1 + 1.0 / npoly) - 2) +
                k * pow(particles[j].density, (1 + 1.0 / npoly) - 2)
            );
            delkern = gradkernel(
                particles[i].x - particles[j].x,
                particles[i].y - particles[j].y,
                h
            );
            particles[i].ax += pressure * delkern;
            particles[j].ax -= pressure * delkern;
            delkern = gradkernel(
                particles[i].y - particles[j].y,
                particles[i].x - particles[j].x,
                h
            );
            particles[i].ay += pressure * delkern;
            particles[j].ay -= pressure * delkern;
        }

        particles[i].vx = particles[i].vx0 + particles[i].ax * dt;
        particles[i].vy = particles[i].vy0 + particles[i].ay * dt;

        particles[i].x += particles[i].vx * dt;
        particles[i].y += particles[i].vy * dt;

        velocity = particles[i].vx;
        particles[i].vx = 0.5 * (particles[i].vx + particles[i].vx0);
        particles[i].vx0 = velocity;
        velocity = particles[i].vy;
        particles[i].vy = 0.5 * (particles[i].vy + particles[i].vy0);
        particles[i].vy0 = velocity;
    }
}

void System::Plot(Canvas* canvas)
{
    canvas->DrawPoints(
        &particles[0],
        particles.size(),
        -1.5,
        1.5,
        -1.5,
        1.5,
        255,
        200,
        50
    );
}

void System::Write(FILE* output)
{
    fprintf(output, "%g, %g, %g, %d, %d, %g, %g, %g\n", lambda, k, m, Npart, npoly, dt, h, damp);
    for (i = 0; i < particles.size(); i++) {
        fprintf(output, "%g, %g, %g\n", particles[i].x, particles[i].y, particles[i].density);
    }
}

int main(int argc, char** argv)
{
    unsigned int WIDTH = 800;
    unsigned int HEIGHT = 600;
    unsigned int N_PART = 500;
    double RADIUS = 1;
    bool running = true;
    clock_t clockCounter;
    FILE* output = NULL;
    Canvas canvas(WIDTH, HEIGHT);
    System system(N_PART, RADIUS);

    if (0 != canvas.GetError()) {
        return 1;
    }

    output = fopen("sphstar.csv", "w");
    if (NULL == output) {
        fprintf(stdout, "Unable to open file.\n");
        return 1;
    }

    while (running) {
        clockCounter = clock();
        running = canvas.HandleEvents();

        while (CLOCKS_PER_SEC / 60.0 > clock() - clockCounter) {
            system.Run();
        }

        canvas.Clear(0, 0, 0);

        system.Plot(&canvas);

        canvas.Show();

        if (!running) {
            system.Write(output);
        }

        fprintf(
            stdout,
            "%.2g\r",
            ((float) CLOCKS_PER_SEC) / (clock() - clockCounter)
        );
    }

    fclose(output);

    return 0;
}
