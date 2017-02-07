#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <vector>
#include "canvas.h"

double powerN(double x, unsigned int n)
{
    double y = 1;
    for (; n > 0; n--) {
        y *= x;
    }
    return y;
}

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
    double nu;
    double damp;
    double kpress;
    double rho0;
    double K;
    double g;
    double boundary;
    bool dam;
    unsigned int Ncells;
    Particle *particles;
    ParticleGrid *grid;
    protected:
        void SetParams();
        void Init();
    public:
        System();
        System(unsigned int);
        ~System();
        void Run();
        void Release();
        void Plot(Canvas*);
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
    SetParams();
    Init();
}

System::~System()
{
    delete particles;
    delete grid;
}

void System::SetParams()
{
    dt = 0.003;
    h = 0.03;
    mass = 10.0 / Npart;
    nu = 0.005;
    damp = 8.0;
    kpress = 0.5;
    rho0 = 3.0;
    K = 5000.0;
    g = -9.8;
    boundary = 2;
    dam = true;
}

void System::Init()
{
    unsigned int totalNcells = 0;

    Ncells = 1.05 * boundary / h;
    totalNcells = Ncells * Ncells;

    particles = new Particle[Npart];
    grid = new ParticleGrid(totalNcells, Npart);

    srand(time(NULL));
    for (i = 0; i < Npart; i++) {
        particles[i].x = (double) rand() / RAND_MAX * boundary / 2 - boundary;
        particles[i].y = (double) rand() / RAND_MAX * 2 * boundary - boundary;
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
    double viscosity = 0;
    double delkernx = 0, delkerny = 0;
    double xij = 0, yij = 0;
    double velocity = 0;
    double xcells = 0, ycells = 0;
    unsigned int totalNcells = Ncells * Ncells;
    Particle *neighbor;

    grid->clear();

    for (i = 0; i < Npart; i++) {
        xcells = (particles[i].x + 1.05 * boundary) / (2.1 * boundary);
        ycells = (particles[i].y + 1.05 * boundary) / (2.1 * boundary);
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

        particles[i].ax = -damp * particles[i].vx;
        if ((dam ? - boundary / 2 : boundary) < particles[i].x) {
            particles[i].ax -= K * (particles[i].x - (dam ? - boundary / 2 : boundary));
        } else if (-boundary > particles[i].x) {
            particles[i].ax -= K * (particles[i].x + boundary);
        }
        particles[i].ay = -damp * particles[i].vy;
        particles[i].ay += g;
        if (boundary < particles[i].y) {
            particles[i].ay -= K * (particles[i].y - boundary);
        } else if (-boundary > particles[i].y) {
            particles[i].ay -= K * (particles[i].y + boundary);
        }

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
                                    (powerN(particles[i].density / rho0, 7) - 1) /
                                    powerN(particles[i].density, 2) +
                                    (powerN(neighbor->density / rho0, 7) - 1) /
                                    powerN(neighbor->density, 2)
                                );
                                xij = particles[i].x - neighbor->x;
                                yij = particles[i].y - neighbor->y;
                                delkernx = gradkernel(xij, yij, h);
                                delkerny = gradkernel(yij, xij, h);
                                viscosity = 2 * nu * mass * (
                                    (xij) * delkernx + (yij) * delkerny
                                ) / (xij*xij + yij*yij + 0.01 * h * h);

                                xij = particles[i].vx - neighbor->vx;
                                yij = particles[i].vy - neighbor->vy;
                                particles[i].ax += pressure * delkernx + (xij) * viscosity / particles[i].density;
                                particles[i].ay += pressure * delkerny + (yij) * viscosity / particles[i].density;
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

void System::Release()
{
    dam = !dam;
    if (dam) {
        damp = 8.0;
    } else {
        damp = 0.05;
    }
}

void System::Plot(Canvas* canvas)
{
    canvas->DrawPoints(
        &particles[0],
        Npart,
        -1.1 * boundary,
        1.1 * boundary,
        -1.1 * boundary,
        1.1 * boundary,
        0,
        0,
        255
    );
}

int main(int argc, char** argv)
{
    unsigned int WIDTH = 800;
    unsigned int HEIGHT = 600;
    unsigned int N_PART = 3000;
    int handle_code = -1;
    bool running = true;
    clock_t clockCounter;
    Canvas canvas(WIDTH, HEIGHT);
    System system(N_PART);

    if (0 != canvas.GetError()) {
        return 1;
    }


    while (running) {
        clockCounter = clock();
        handle_code = canvas.HandleEvents();
        if (0 == handle_code) {
            running = false;
        } else if (1 == handle_code) {
            system.Release();
        }

        while (CLOCKS_PER_SEC / 60.0 > clock() - clockCounter) {
            system.Run();
        }

        canvas.Clear();

        system.Plot(&canvas);

        canvas.Show();

        fprintf(
            stdout,
            "%.2g \r",
            ((float) CLOCKS_PER_SEC) / (clock() - clockCounter)
        );
        fflush(stdout);
    }

    return 0;
}
