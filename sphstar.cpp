#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <vector>
#include "canvas.h"

using namespace std;

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

typedef struct Particle
{
    double x, y;
    double vx, vy;
    double vx0, vy0;
    double ax, ay;
    double density;
    unsigned int cell;
}Particle;

class ParticleGrid
{
    unsigned int i, j;
    unsigned int grid_size;
    unsigned int total_size;
    unsigned int *actual_size;
    Particle ***particles;
    public:
        ParticleGrid(unsigned int, unsigned int);
        ~ParticleGrid();
        Particle**& operator[](unsigned int);
        unsigned int size(unsigned int);
        int push_back(unsigned int, Particle*);
        void clear();
};

ParticleGrid::ParticleGrid(unsigned int grid_size_to_allocate, unsigned int size_to_allocate)
{
    grid_size = grid_size_to_allocate;
    total_size = size_to_allocate;
    particles = new Particle**[grid_size];
    actual_size = new unsigned int[grid_size];
    for (i = 0; i < grid_size; i++) {
        particles[i] = new Particle*[total_size];
        for (j = 0; j < total_size; j++) {
            particles[i][j] = NULL;
        }
        actual_size[i] = 0;
    }
}

ParticleGrid::~ParticleGrid()
{
    for (i = 0; i < grid_size; i++) {
        delete [] particles[i];
    }
    delete [] particles;
}

Particle**& ParticleGrid::operator[](unsigned int index)
{
    return particles[index];
}

unsigned int ParticleGrid::size(unsigned int index)
{
    return actual_size[index];
}

int ParticleGrid::push_back(unsigned int index, Particle* particle)
{
    if (actual_size[index] < total_size) {
        particles[index][actual_size[index]] = particle;
        actual_size[index]++;
        return 0;
    }
    return 1;
}

void ParticleGrid::clear()
{
    for (i = 0; i < grid_size; i++) {
        /*for (j = 0; j < actual_size[i]; j++) {
            particles[i][j] = NULL;
        }*/
        actual_size[i] = 0;
    }
}

class System
{
    unsigned int Npart;
    unsigned int i, j;
    int m, n;
    double dt;
    double h;
    double mass;
    double damp;
    double kpress;
    double R0;
    int npoly;
    double lambda;
    double boundary;
    unsigned int Ncells;
    Particle *particles;
    ParticleGrid *grid;
    protected:
        void SetParams();
        void Init();
    public:
        System();
        System(unsigned int, double);
        ~System();
        void Run();
        void Plot(Canvas*);
        void Write(FILE*);
};

System::System()
{
    Npart = 100;
    R0 = 5;
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
    delete grid;
    delete particles;
}

void System::SetParams()
{
    dt = 0.01;
    h = 0.1;
    mass = 1.33333 / Npart;
    damp = 1.0;
    kpress = 0.25;
    npoly = 1;
    lambda = 1.0;
    boundary = 2.0;
}

void System::Init()
{
    double radius = 0;
    double phi = 0;
    unsigned int totalNcells;

    Ncells = boundary / h;
    totalNcells = Ncells * Ncells;

    grid = new ParticleGrid(totalNcells, Npart);

    particles = new Particle[Npart];

    srand(time(NULL));
    for (i = 0; i < Npart; i++) {
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
        particles[i].cell = 0;
    }
}

void System::Run()
{
    double rho = 0;
    double pressure = 0;
    double delkern = 0;
    double velocity = 0;
    double xcells, ycells;
    unsigned int totalNcells = Ncells * Ncells;
    Particle *neighbor;

    grid->clear();

    for (i = 0; i < Npart; i++) {
        xcells = 0.5 * (particles[i].x + boundary) / boundary;
        ycells = 0.5 * (particles[i].y + boundary) / boundary;
        if (
            1.0 > xcells && 0 <= xcells &&
            1.0 > ycells && 0 <= ycells
        ) {
            particles[i].cell = Ncells * (
                xcells +
                (int) (ycells * Ncells)
            );
        } else {
            particles[i].cell = totalNcells;
        }
        if (totalNcells > particles[i].cell) {
            grid->push_back(particles[i].cell, &particles[i]);
        }
    }

    for (i = 0; i < Npart; i++) {
        particles[i].density = mass * kernel(0, 0, h);
        if (totalNcells > particles[i].cell) {
            for (m = -1; m < 2; m++) {
                for (n = -1; n < 2; n++) {
                    if (
                        !(0 == particles[i].cell / Ncells && -1 == m) &&
                        !(0 == particles[i].cell % Ncells && -1 == n) &&
                        !(Ncells - 1 == particles[i].cell / Ncells && 1 == m) &&
                        !(Ncells - 1 == particles[i].cell % Ncells && 1 == n)

                    ) {
                        for (j = 0; j < grid->size(particles[i].cell + n + Ncells * m); j++) {
                            neighbor = (*grid)[particles[i].cell + n + Ncells * m][j];
                            if (NULL != neighbor) {
                                rho = mass * kernel(
                                    particles[i].x - neighbor->x,
                                    particles[i].y - neighbor->y,
                                    h
                                );
                                particles[i].density += rho;
                            }
                        }
                    }
                }
            }
        }
    }

    for (i = 0; i < Npart; i++) {
        particles[i].ax = - damp * particles[i].vx - lambda * particles[i].x;
        particles[i].ay = - damp * particles[i].vy - lambda * particles[i].y;

        if (totalNcells > particles[i].cell) {
            for (m = -1; m < 2; m++) {
                for (n = -1; n < 2; n++) {
                    if (
                        !(0 == particles[i].cell / Ncells && -1 == m) &&
                        !(0 == particles[i].cell % Ncells && -1 == n) &&
                        !(Ncells - 1 == particles[i].cell / Ncells && 1 == m) &&
                        !(Ncells - 1 == particles[i].cell % Ncells && 1 == n)

                    ) {
                        for (j = 0; j < grid->size(particles[i].cell + n + Ncells * m); j++) {
                            neighbor = (*grid)[particles[i].cell + n + Ncells * m][j];
                            if (NULL != neighbor) {
                                pressure = -mass * kpress * (
                                    pow(particles[i].density, (1 + 1.0 / npoly) - 2) +
                                    pow(neighbor->density, (1 + 1.0 / npoly) - 2)
                                );
                                delkern = gradkernel(
                                    particles[i].x - neighbor->x,
                                    particles[i].y - neighbor->y,
                                    h
                                );
                                particles[i].ax += pressure * delkern;
                                delkern = gradkernel(
                                    particles[i].y - neighbor->y,
                                    particles[i].x - neighbor->x,
                                    h
                                );
                                particles[i].ay += pressure * delkern;
                            }
                        }
                    }
                }
            }
        }
    }

    for (i = 0; i < Npart; i++) {
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
        Npart,
        -boundary,
        boundary,
        -boundary,
        boundary,
        255,
        200,
        50
    );
}

void System::Write(FILE* output)
{
    fprintf(output, "%g, %g, %g, %d, %d, %g, %g, %g\n", lambda, kpress, mass, Npart, npoly, dt, h, damp);
    for (i = 0; i < Npart; i++) {
        fprintf(output, "%g, %g, %g\n", particles[i].x, particles[i].y, particles[i].density);
    }
}

int main(int argc, char** argv)
{
    unsigned int WIDTH = 800;
    unsigned int HEIGHT = 600;
    unsigned int N_PART = 1000;
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
            "%.2g \r",
            ((float) CLOCKS_PER_SEC) / (clock() - clockCounter)
        );
        fflush(stdout);
    }

    fclose(output);

    return 0;
}
