#include "canvas.h"

using namespace std;

Canvas::Canvas()
{
    xScale = 1;
    yScale = 1;
    xTrans = 0;
    yTrans = 0;
    width = 640;
    height = 480;
    Init();
}

Canvas::Canvas(unsigned int w, unsigned int h)
{
    xScale = 1;
    yScale = 1;
    xTrans = 0;
    yTrans = 0;
    width = w;
    height = h;
    Init();
}

Canvas::~Canvas()
{
    Quit();
}

void Canvas::Init()
{
    if (0 > SDL_Init(SDL_INIT_VIDEO)) {
        fprintf(stderr,"SDL_Init failed:\n%s\n", SDL_GetError());
        error_code = 1;
    }

    display = SDL_CreateWindow(
        "SPH",
        SDL_WINDOWPOS_UNDEFINED,
        SDL_WINDOWPOS_UNDEFINED,
        width,
        height,
        0
    );
    if (NULL == display) {
        fprintf(stderr, "SDL_CreateWindow failed:\n%s\n", SDL_GetError());
        error_code = 2;
    }

    renderer = SDL_CreateRenderer(display, -1, 0);
    if (NULL == renderer) {
        fprintf(stderr, "SDL_CreateRenderer failed:\n%s\n", SDL_GetError());
        error_code = 3;
    }

    error_code = 0;
}

unsigned short int Canvas::HandleEvents()
{
    while (SDL_PollEvent(&event)) {
        if (SDL_QUIT == event.type) {
            return 0;
        } else if (SDL_KEYUP == event.type) {
            if (SDLK_ESCAPE == event.key.keysym.sym) {
                return 0;
            } else if (SDLK_SPACE == event.key.keysym.sym) {
                return 1;
            }
        }
    }
    return -1;
}

void Canvas::Clear(unsigned short r, unsigned short g, unsigned short b)
{
    SDL_SetRenderDrawColor(renderer, r, g, b, 255);
    SDL_RenderClear(renderer);
}

void Canvas::Show()
{
    SDL_RenderPresent(renderer);
}

void Canvas::Quit()
{
    SDL_DestroyWindow(display);
    SDL_Quit();
}

int Canvas::GetError()
{
    return error_code;
}

void Canvas::CalculateTransformParam(
    double leftX,
    double rightX,
    double bottomY,
    double topY
) {
    yScale =  height / (bottomY - topY);
    xScale =  rightX > leftX ? fabs(yScale) : -fabs(yScale);
    xTrans =  width / 2.0 - xScale * (leftX + rightX) / 2;
    yTrans =  height / 2.0 - yScale * (bottomY + topY) / 2;
}
