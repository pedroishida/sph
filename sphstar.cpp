#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <vector>
#include <thread>
#include "canvas.h"

using namespace std;

double kernel(double x, double y, double z, double h)
{
    double q = sqrt(x*x + y*y + z*z) / h;
    double C = 1.0 / (4 * M_PI * h*h*h);
    if (0 <= q && 1 > q) {
        return C * (4 - 6 * q*q + 3 * q*q*q);
    } else if (1 <= q && 2 > q) {
        return C * (8 - 12 * q + 6 * q*q - q*q*q);
    } else {
        return 0;
    }
}

double gradkernel(double x, double y, double z, double h)
{
    double q = sqrt(x*x + y*y + z*z) / h;
    double C = 1.0 / (4 * M_PI * h*h*h);
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
    double x, y, z;
    double vx, vy, vz;
    double vx0, vy0, vz0;
    double ax, ay, az;
    double density;
    unsigned int cell;
}Particle;

class ParticleGrid
{
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
    for (unsigned int i = 0; i < grid_size; i++) {
        particles[i] = new Particle*[total_size];
        for (unsigned int j = 0; j < total_size; j++) {
            particles[i][j] = NULL;
        }
        actual_size[i] = 0;
    }
}

ParticleGrid::~ParticleGrid()
{
    for (unsigned int i = 0; i < grid_size; i++) {
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
    for (unsigned int i = 0; i < grid_size; i++) {
        /*for (unsigned int j = 0; j < actual_size[i]; j++) {
            particles[i][j] = NULL;
        }*/
        actual_size[i] = 0;
    }
}

class System
{
    unsigned int Nthreads;
    unsigned int Npart;
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
    thread *calcThreads;
    Particle *particles;
    ParticleGrid *grid;
    protected:
        void SetParams();
        void Init();
        void Randomize(unsigned int, unsigned int);
        void CalculateCell(unsigned int, unsigned int);
        void CalculateDensity(unsigned int, unsigned int);
        void CalculateAcceleration(unsigned int, unsigned int);
        void CalculatePosition(unsigned int, unsigned int);
    public:
        System();
        System(unsigned int, double, unsigned int);
        ~System();
        int Run();
        void Plot(Canvas*);
        void Write(FILE*);
};

System::System()
{
    Nthreads = 1;
    Npart = 100;
    R0 = 5;
    SetParams();
    Init();
}

System::System(unsigned int np, double radius0, unsigned int nt)
{
    Nthreads = nt;
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
    h = 0.15;
    mass = 1.33333 / Npart;
    damp = 1.0;
    kpress = 0.25;
    npoly = 1;
    lambda = 1.0;
    boundary = 2.0;
}

void System::Init()
{
    unsigned int totalNcells;
    unsigned int dNpart = Npart / Nthreads;

    calcThreads = new thread[Nthreads];

    Ncells = boundary / h;
    totalNcells = Ncells * Ncells * Ncells;

    grid = new ParticleGrid(totalNcells, Npart);

    particles = new Particle[Npart];

    srand(time(NULL));

    for (unsigned int i = 0; i < Nthreads; i++) {
        calcThreads[i] = thread(&System::Randomize, this, i * dNpart, i == Nthreads - 1 ? Npart : (i + 1) * dNpart);
    }
    for (unsigned int i = 0; i < Nthreads; i++) {
        calcThreads[i].join();
    }

}

void System::Randomize(unsigned int start, unsigned int end)
{
    double radius = 0;
    double costheta = 0;
    double phi = 0;

    for (unsigned int i = start; i < end; i++) {
        radius = pow((double) rand() / RAND_MAX, 1.0 / 3.0) * R0;
        costheta = ((double) rand() / RAND_MAX) * 2.0 - 1.0;
        phi = (double) rand() / RAND_MAX * 2 * M_PI;
        particles[i].x = cos(phi) * sqrt(1 - costheta*costheta) * radius;
        particles[i].y = sin(phi) * sqrt(1 - costheta*costheta) * radius;
        particles[i].z = costheta * radius;
        particles[i].vx = 0;
        particles[i].vy = 0;
        particles[i].vz = 0;
        particles[i].vx0 = 0;
        particles[i].vy0 = 0;
        particles[i].vz0 = 0;
        particles[i].density = 0;
        particles[i].ax = 0;
        particles[i].ay = 0;
        particles[i].az = 0;
        particles[i].cell = 0;
    }
}

int System::Run()
{
    unsigned int dNpart = Npart / Nthreads;

    grid->clear();

    for (unsigned int i = 0; i < Nthreads; i++) {
        calcThreads[i] = thread(&System::CalculateCell, this, i * dNpart, i == Nthreads - 1 ? Npart : (i + 1) * dNpart);
    }
    for (unsigned int i = 0; i < Nthreads; i++) {
        calcThreads[i].join();
    }

    for (unsigned int i = 0; i < Nthreads; i++) {
        calcThreads[i] = thread(&System::CalculateDensity, this, i * dNpart, i == Nthreads - 1 ? Npart : (i + 1) * dNpart);
    }
    for (unsigned int i = 0; i < Nthreads; i++) {
        calcThreads[i].join();
    }

    for (unsigned int i = 0; i < Nthreads; i++) {
        calcThreads[i] = thread(&System::CalculateAcceleration, this, i * dNpart, i == Nthreads - 1 ? Npart : (i + 1) * dNpart);
    }
    for (unsigned int i = 0; i < Nthreads; i++) {
        calcThreads[i].join();
    }

    for (unsigned int i = 0; i < Nthreads; i++) {
        calcThreads[i] = thread(&System::CalculatePosition, this, i * dNpart, i == Nthreads - 1 ? Npart : (i + 1) * dNpart);
    }
    for (unsigned int i = 0; i < Nthreads; i++) {
        calcThreads[i].join();
    }

    return 0;
}

void System::CalculateCell(unsigned int start, unsigned int end)
{
    double xcells, ycells, zcells;
    unsigned int totalNcells = Ncells * Ncells * Ncells;

    for (unsigned int i = start; i < end; i++) {
        xcells = 0.5 * (particles[i].x + boundary) / boundary;
        ycells = 0.5 * (particles[i].y + boundary) / boundary;
        zcells = 0.5 * (particles[i].z + boundary) / boundary;
        if (
            1.0 > xcells && 0 <= xcells &&
            1.0 > ycells && 0 <= ycells &&
            1.0 > zcells && 0 <= zcells
        ) {
            particles[i].cell = Ncells * (
                xcells +
                (int) (ycells * Ncells) +
                ((int) (zcells * Ncells)) * Ncells
            );
        } else {
            particles[i].cell = totalNcells;
        }
        if (totalNcells > particles[i].cell) {
            grid->push_back(particles[i].cell, &particles[i]);
        }
    }
}

void System::CalculateDensity(unsigned int start, unsigned int end)
{
    double rho = 0;
    unsigned int totalNcells = Ncells * Ncells * Ncells;
    unsigned int partialNcells = Ncells * Ncells;
    Particle *neighbor;

    for (unsigned int i = start; i < end; i++) {
        particles[i].density = mass * kernel(0, 0, 0, h);
        if (totalNcells > particles[i].cell) {
            for (int l = -1; l < 2; l++) {
                for (int m = -1; m < 2; m++) {
                    for (int n = -1; n < 2; n++) {
                        if (
                            !(0 == (particles[i].cell / (partialNcells)) && -1 == l) &&
                            !(0 == (particles[i].cell % (partialNcells)) / Ncells && -1 == m) &&
                            !(0 == (particles[i].cell % (partialNcells)) % Ncells && -1 == n) &&
                            !(Ncells - 1 == (particles[i].cell / (partialNcells)) && 1 == l) &&
                            !(Ncells - 1 == (particles[i].cell % (partialNcells)) / Ncells && 1 == m) &&
                            !(Ncells - 1 == (particles[i].cell % (partialNcells)) % Ncells && 1 == n)

                        ) {
                            for (unsigned int j = 0; j < grid->size(particles[i].cell + n + Ncells * m + partialNcells * l); j++) {
                                neighbor = (*grid)[particles[i].cell + n + Ncells * m + partialNcells * l][j];
                                if (NULL != neighbor) {
                                    rho = mass * kernel(
                                        particles[i].x - neighbor->x,
                                        particles[i].y - neighbor->y,
                                        particles[i].z - neighbor->z,
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
    }
}

void System::CalculateAcceleration(unsigned int start, unsigned int end)
{
    double delkern = 0;
    double pressure = 0;
    unsigned int totalNcells = Ncells * Ncells * Ncells;
    unsigned int partialNcells = Ncells * Ncells;
    Particle *neighbor;

    for (unsigned int i = start; i < end; i++) {
        particles[i].ax = - damp * particles[i].vx - lambda * particles[i].x;
        particles[i].ay = - damp * particles[i].vy - lambda * particles[i].y;
        particles[i].az = - damp * particles[i].vz - lambda * particles[i].z;

        if (totalNcells > particles[i].cell) {
            for (int l = -1; l < 2; l++) {
                for (int m = -1; m < 2; m++) {
                    for (int n = -1; n < 2; n++) {
                        if (
                            !(0 == (particles[i].cell / (partialNcells)) && -1 == l) &&
                            !(0 == (particles[i].cell % (partialNcells)) / Ncells && -1 == m) &&
                            !(0 == (particles[i].cell % (partialNcells)) % Ncells && -1 == n) &&
                            !(Ncells - 1 == (particles[i].cell / (partialNcells)) && 1 == l) &&
                            !(Ncells - 1 == (particles[i].cell % (partialNcells)) / Ncells && 1 == m) &&
                            !(Ncells - 1 == (particles[i].cell % (partialNcells)) % Ncells && 1 == n)

                        ) {
                            for (unsigned int j = 0; j < grid->size(particles[i].cell + n + Ncells * m + partialNcells * l); j++) {
                                neighbor = (*grid)[particles[i].cell + n + Ncells * m + partialNcells * l][j];
                                if (NULL != neighbor) {
                                    pressure = -mass * kpress * (
                                        pow(particles[i].density, (1 + 1.0 / npoly) - 2) +
                                        pow(neighbor->density, (1 + 1.0 / npoly) - 2)
                                    );
                                    delkern = gradkernel(
                                        particles[i].x - neighbor->x,
                                        particles[i].y - neighbor->y,
                                        particles[i].z - neighbor->z,
                                        h
                                    );
                                    particles[i].ax += pressure * delkern;
                                    delkern = gradkernel(
                                        particles[i].y - neighbor->y,
                                        particles[i].x - neighbor->x,
                                        particles[i].z - neighbor->z,
                                        h
                                    );
                                    particles[i].ay += pressure * delkern;
                                    delkern = gradkernel(
                                        particles[i].z - neighbor->z,
                                        particles[i].y - neighbor->y,
                                        particles[i].x - neighbor->x,
                                        h
                                    );
                                    particles[i].az += pressure * delkern;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void System::CalculatePosition(unsigned int start, unsigned int end)
{
    double velocity;

    for (unsigned int i = start; i < end; i++) {
        particles[i].vx = particles[i].vx0 + particles[i].ax * dt;
        particles[i].vy = particles[i].vy0 + particles[i].ay * dt;
        particles[i].vz = particles[i].vz0 + particles[i].az * dt;

        particles[i].x += particles[i].vx * dt;
        particles[i].y += particles[i].vy * dt;
        particles[i].z += particles[i].vz * dt;

        velocity = particles[i].vx;
        particles[i].vx = 0.5 * (particles[i].vx + particles[i].vx0);
        particles[i].vx0 = velocity;
        velocity = particles[i].vy;
        particles[i].vy = 0.5 * (particles[i].vy + particles[i].vy0);
        particles[i].vy0 = velocity;
        velocity = particles[i].vz;
        particles[i].vz = 0.5 * (particles[i].vz + particles[i].vz0);
        particles[i].vz0 = velocity;
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
    fprintf(output, "lambda,kpress,mass,Npart,npoly,dt,h,damp\n");
    fprintf(output, "%g,%g,%g,%d,%d,%g,%g,%g\n", lambda, kpress, mass, Npart, npoly, dt, h, damp);
    fprintf(output, "x,y,z,density\n");
    for (unsigned int i = 0; i < Npart; i++) {
        fprintf(output, "%g,%g,%g,%g\n", particles[i].x, particles[i].y, particles[i].z, particles[i].density);
    }
}

int main(int argc, char** argv)
{
    unsigned int THREADS = 4;
    unsigned int WIDTH = 800;
    unsigned int HEIGHT = 600;
    unsigned int N_PART = 5000;
    double RADIUS = 1;
    bool running = true;
    int stable = 0;
    clock_t clockCounter;
    FILE* output = NULL;
    Canvas canvas(WIDTH, HEIGHT);
    System system(N_PART, RADIUS, THREADS);

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
            stable = system.Run();
        }

        canvas.Clear(0, 0, 0);

        system.Plot(&canvas);

        canvas.Show();

        if (!running) {
            system.Write(output);
        }

        if (stable) {
            fprintf(stdout, "STABLE ");
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
