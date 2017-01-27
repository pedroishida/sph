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
    vector<SDL_Point> points;
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
    CalculateTransformParam(leftX, rightX, bottomY, topY);

    if (0 < size && NULL != positions) {
        points.resize(5 * size);

        for (i = 0; i < points.size(); i+=5) {
            j = i / 5;
            points[i].x = positions[j].x * xScale + xTrans;
            points[i].y = positions[j].y * yScale + yTrans;
            points[i + 1].x = points[i].x + 1;
            points[i + 1].y = points[i].y;
            points[i + 2].x = points[i].x - 1;
            points[i + 2].y = points[i].y;
            points[i + 3].x = points[i].x;
            points[i + 3].y = points[i].y + 1;
            points[i + 4].x = points[i].x;
            points[i + 4].y = points[i].y - 1;
        }

        SDL_SetRenderDrawColor(renderer, r, g, b, 255);
        SDL_RenderDrawPoints(renderer, &points[0], points.size());
    }
}

#endif
