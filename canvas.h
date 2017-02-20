#ifndef CANVAS_H
#define CANVAS_H

#include <cstdio>
#include <cmath>
#include <vector>
#include <SDL2/SDL.h>

using namespace std;

class Canvas
{
    unsigned int i, j;
    unsigned int width;
    unsigned int height;
    double xScale, yScale;
    double xTrans, yTrans;
    int error_code;
    SDL_Window* display;
    SDL_Renderer* renderer;
    SDL_Event event;
    protected:
        void Init();
        void Quit();
    public:
        Canvas();
        Canvas(unsigned int, unsigned int);
        ~Canvas();
        unsigned short int HandleEvents();
        void Clear(
            unsigned short int = 255,
            unsigned short int = 255,
            unsigned short int = 255
        );
        void Show();
        int GetError();
        void CalculateTransformParam(double, double, double, double);
        template <typename xy>
        void DrawPoints(
            xy*,
            unsigned int,
            double,
            double,
            double,
            double,
            unsigned short int = 0,
            unsigned short int = 0,
            unsigned short int = 0
        );
};

template <typename xy>
void Canvas::DrawPoints(
    xy* positions,
    unsigned int size,
    double leftX,
    double rightX,
    double bottomY,
    double topY,
    unsigned short int r,
    unsigned short int g,
    unsigned short int b
) {
    double x, y;

    if (0 < size && NULL != positions) {
        CalculateTransformParam(leftX, rightX, bottomY, topY);
        SDL_SetRenderDrawColor(renderer, r, g, b, 255);
        for (i = 0; i < size; i++) {
            x = positions[i].x * xScale + xTrans;
            y = positions[i].y * yScale + yTrans;
            SDL_RenderDrawPoint(renderer, x, y);
            SDL_RenderDrawPoint(renderer, x + 1, y);
            SDL_RenderDrawPoint(renderer, x - 1, y);
            SDL_RenderDrawPoint(renderer, x, y + 1);
            SDL_RenderDrawPoint(renderer, x, y - 1);
        }
    }
}

#endif
