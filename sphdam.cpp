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

class System
{
    unsigned int Npart;
    unsigned int i, j;
    double dt;
    double h;
    double m;
    double nu;
    double damp;
    double kpress;
    double rho0;
    double K;
    double g;
    double box;
    bool dam;
    vector<Particle> particles;
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
}

void System::SetParams()
{
    dt = 0.001666;
    h = 0.0666;
    m = 5.0 / Npart;
    nu = 0.005;
    damp = 8.0;
    kpress = 0.5;
    rho0 = 2.861;
    K = 5000.0;
    g = -9.8;
    box = 1.5;
    dam = true;
}

void System::Init()
{
    particles.resize(Npart);

    srand(time(NULL));
    for (i = 0; i < particles.size(); i++) {
        particles[i].x = (double) rand() / RAND_MAX * box / 2 - box;
        particles[i].y = (double) rand() / RAND_MAX * -box;
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

        particles[i].ax -= damp * particles[i].vx;
        if (dam) {
            if (-box / 2 < particles[i].x) {
                particles[i].ax -= K * (particles[i].x + box / 2);
            }
        } else {
            if (box < particles[i].x) {
                particles[i].ax -= K * (particles[i].x - box);
            }
        }
        if (-box > particles[i].x) {
            particles[i].ax -= K * (particles[i].x + box);
        }
        particles[i].ay -= damp * particles[i].vy;
        particles[i].ay += g;
        if (box < particles[i].y) {
            particles[i].ay -= K * (particles[i].y - box);
        } else if (-box > particles[i].y) {
            particles[i].ay -= K * (particles[i].y + box);
        }

        for (j = i + 1; j < particles.size(); j++) {
            pressure = -m * kpress * (
                (powerN(particles[i].density / rho0, 7) - 1) /
                (particles[i].density*particles[i].density) +
                (powerN(particles[j].density / rho0, 7) - 1) /
                (particles[j].density*particles[j].density)
            );
            xij = particles[i].x - particles[j].x;
            yij = particles[i].y - particles[j].y;
            delkernx = gradkernel(xij, yij, h);
            delkerny = gradkernel(yij, xij, h);
            viscosity = 2 * nu * m * (
                (xij) * delkernx + (yij) * delkerny
            ) / (powerN(xij, 2) + powerN(yij, 2) + 0.01 * h * h);

            xij = particles[i].vx - particles[j].vx;
            yij = particles[i].vy - particles[j].vy;
            particles[i].ax += pressure * delkernx + (xij) * viscosity / particles[j].density;
            particles[j].ax -= pressure * delkernx + (xij) * viscosity / particles[j].density;
            particles[i].ay += pressure * delkerny + (yij) * viscosity / particles[i].density;
            particles[j].ay -= pressure * delkerny + (yij) * viscosity / particles[i].density;
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
        particles.size(),
        -1.5 * box,
        1.5 * box,
        -1.5 * box,
        1.5 * box,
        0,
        0,
        255
    );
}

int main(int argc, char** argv)
{
    unsigned int WIDTH = 800;
    unsigned int HEIGHT = 600;
    unsigned int N_PART = 300;
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
            "%.2g\r",
            ((float) CLOCKS_PER_SEC) / (clock() - clockCounter)
        );
    }

    return 0;
}
